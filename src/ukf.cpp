#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

    /*
     * if this is false, laser measurements will be ignored (except during
     * init)
     */
    use_laser_ = true;

    /*
     * if this is false, radar measurements will be ignored (except during 
     * init)
     */
    use_radar_ = true;


    /*
     * initial state vector, 5 state in CTRV model.
     * px, py, v, yaw, yaw rate
     */
    x_ = VectorXd(5);

    /*
     * initial covariance matrix, 5x5
     */
    P_ = MatrixXd(5, 5);

    /*
     * Process noise standard deviation longitudinal acceleration in m/s^2
     *
     * Assuming Max acceleration noise parameter for bicycle with a smooth
     * response is 1 m/sec and taking 1/2 of it as acceleration noise
     * standard deviation here.
     */
    std_a_ = 0.5 ;

    /*
     * Process noise standard deviation yaw acceleration in rad/s^2
     * Assuming an yaw acceleration of phi/6 ~= 0.5 rad/sec^2
     * So for a bicycle traveling in a circle with a const yaw rate of
     * let say phi/8 rad/sec which will take (phi/8)*16s = 2phi
     * => 16 second to complete the circle.
     * Now Assume an yaw acceleration of -phi/6 so in 1 second the
     * radial velocity will be phi/8 - phi/6*1 = -phi/24. To complete
     * the full circle in other direction will now take 48 second which
     * "I think" is a valid smooth transition in other direction and
     * is assumed here.
     *
     */
    std_yawdd_ = 0.5;

    /*
     * Laser measurement noise standard deviation position1 in m
     */
    std_laspx_ = 0.15;

    /*
     * Laser measurement noise standard deviation position2 in m
     */
    std_laspy_ = 0.15;

    /*
     * Radar measurement noise standard deviation radius in m
     */
    std_radr_ = 0.3;

    /*
     * Radar measurement noise standard deviation angle in rad
     */
    std_radphi_ = 0.03;

    /*
     * Radar measurement noise standard deviation radius change in m/s
     */
    std_radrd_ = 0.3;

    /*
     * state dimention in CTRV (Const rate and velocity model)
     */
    n_x_ = 5;

    /*
     * 5 state + 2 process noise (longitudinal acceleration and
     * yow acceleration)
     */
    n_aug_ = 7;

    /*
     * spreading parameter
     */
    lambda_ = 3 - n_aug_;

    /*
     * Predicted sigma point matrix
     */
    Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

    /*
     * Vector for weights
     */
    weights_ = VectorXd(2*n_aug_ + 1);

    /*
     * Augmented mean vector
     */
     x_arg = VectorXd(n_aug_);

    /*
     * Augmented state covariance matrix
     */
    P_arg = MatrixXd(n_aug_, n_aug_);

    /*
     * Sigma point matrix
     */
    Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    /*
     * Process noise covariance matrix
     */
    Q = MatrixXd(2, 2);

    Q << std_a_*std_a_, 0,
         0, std_yawdd_*std_yawdd_;

    previous_timestamp_ = 0;

    /*
     * Measurement Noise covariance for radar
     * Constant for radar, no need to calculate
     * on every measurement
     */
    R_radar = MatrixXd(3,3);
    R_radar << std_radr_*std_radr_, 0, 0,
                0, std_radphi_*std_radphi_, 0,
                0,  0, std_radrd_*std_radrd_;

    /*
     * Measurement Noise covariance for laser
     * no need to calculate on every measurement
     */
    R_laser = MatrixXd(2,2);
    R_laser << std_laspx_*std_laspx_, 0,
               0, std_laspy_*std_laspy_;


    /*
     * NIS
     */
    NIS_laser_ = 0;
    NIS_radar_ = 0;
}


UKF::~UKF() {}

static double total_time = 0;

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

    int i; 
    float dt; 
    float ro, phi, ro_dot; /* polar coordinate */
    float px, py, vx, vy;  /* cartesian coordinate */

    if (!is_initialized_) {
    
        cout << "Initializing UKF:" <<endl;
        previous_timestamp_ = meas_package.timestamp_;
        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          
            /* first measurement */

            /*
             * Radar data isin polar coordinate
             * Convert to cartesian coordinates and initialize state.
             *
             * Map initial motion state from polar to cartesian  coordinate
             */
            ro = meas_package.raw_measurements_[0];
            phi = meas_package.raw_measurements_[1];
            ro_dot = meas_package.raw_measurements_[2];

            px = ro*cos(phi);
            py = ro*sin(phi);
            
            vx = ro_dot*cos(phi);
            vy = ro_dot*sin(phi);
            x_ << px, py, 5.0, 0, 0;

            /*
             * State covariance matrix Radar
             * - Radar measurement noise standard deviation for measured
             *   radius is 0.3, derive the deviation in x and y direction
             *   as bellow
             */
            double un_x = 0.3*cos(phi);
            double un_y = 0.3*sin(phi);

            /*
             * Radar measurement noise standard deviation angle in rad
             * Given from manufacturing handbook.
             */
            double un_phi = 0.03; 

            /*
             * - initial variance or uncertainty in px and py would assume
             *   within un_x*un_x for px and un_y*un_y for py and
             *   un_phi*un_phi for yaw. Remaining keep upto 1
             */
            P_ << un_x*un_x, 0, 0, 0, 0,
                  0, un_y*un_y, 0, 0, 0,
                  0, 0, 1, 0, 0,
                  0, 0, 0, un_phi*un_phi, 0,
                  0, 0, 0, 0, 1;

            /*
             * if in fused or radar only mode then
             * set the init flag
             */
            if (use_radar_) {
                is_initialized_ = true;
            }
            return;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

           /*
             * Lidar is in cartesian axis and can not measure velocity,
             * so just put vx, vy = 0, 0
             * Initialize motion state.
             */
            px  = meas_package.raw_measurements_[0];
            py = meas_package.raw_measurements_[1]; 
            x_ << px, py, 0, 0, 0;
            
            /*
             * State covariance matrix Lider
             *  - Laser measurement noise standard deviation for px and py is
             *    0.15 given based on the manufacturing handbook
             *  - initial variance or uncertainty in px and py would assume
             *    with in 0.15x0.15 = 0.0225
             *  - keep the uncertainity for v, yaw and yaw rate to 1, higher
             */
            P_ << 0.0225, 0, 0, 0, 0,
                  0, 0.0225, 0, 0, 0,
                  0, 0, 1, 0, 0,
                  0, 0, 0, 1, 0,
                  0, 0, 0, 0, 1;

            /*
             * if in fused or laser only mode then
             * set the init flag
             */
            if (use_laser_) {
                is_initialized_ = true;
            }
            return;
        }
    }

    /*
     * compute the time elapsed between the current and previous measurements
     * dt - expressed in seconds
     */
    dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    if (dt <= 0.001) {
        cout << "dt is very small" << dt << endl;
        return;
    }

    total_time += dt;

    /*
     * Unscented Kalman filter predict
     */
    Prediction(dt);
        
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
        if (use_radar_ == false) {
            return;
        }
        UpdateRadar(meas_package);

        /*
         * Following output is to capture and plot NIS/yaw/yaw rate
         * with time k for radar, to make it simple I just print it
         * to stdout and copy/paste. Data could have been saved to a
         * file  directly but then the responsibility will come to 
         * close the file at exist and currently we closed this program
         * by ctrl + c, at which point should the file be closed. To do
         * that a signal handler can be installed. I just wanted to avoid
         * all these extra noise and keep the program simple focus on
         * solving the ukf. :-)
         */
        cout <<'R'<<'\t'<<NIS_radar_<<'\t'<<total_time<<'\t'<<x_[3]<<'\t'<<yaw_gt;
        cout <<'\t'<<x_[4]<<'\t'<<yawd_gt<<endl;
    } else if (meas_package.sensor_type_ == 
                                        MeasurementPackage::LASER) {
        if (use_laser_ == false) {
            return;
        }
        UpdateLidar(meas_package);
        cout <<'L'<<'\t'<<NIS_laser_<<'\t'<<total_time<<'\t'<<x_[3]<<'\t'<<yaw_gt;
        cout <<'\t'<<x_[4]<<'\t'<<yawd_gt<<endl;
    }
}

/*
 * Prediction
 *
 * Predicts sigma points, the state, and the state covariance matrix.
 * @in {double} delta_t: the change in time (in seconds) between the last
 *                       measurement and this one.
 * @return void.
 */
void UKF::Prediction(double delta_t) {
    int i;
    MatrixXd A;         // square root matrix

    /*
     * Augmented states
     */
    double p_x, p_y, v, yaw, yawd, nu_a, nu_yawdd;

    /*
     * predicted state values
     */
    double px_p, py_p, v_p, yaw_p, yawd_p;


    /*
     * Complete this function! Estimate the object's location. Modify the state
     * vector, x_. Predict sigma points, the state, and the state covariance
     * matrix.
     */

    /********************************************************
     * 0. Generating the sigma points
     *
     ********************************************************/

    /*
     * create augmented mean state
     */
    x_arg.head(5) = x_;
    x_arg(5) = 0;
    x_arg(6) = 0;
    
    /*
     * create augmented covariance matrix
     * Q = process noise covariance matrix
     */
    P_arg.fill(0);
    P_arg.topLeftCorner(5,5) = P_;
    P_arg.bottomRightCorner(2,2) = Q;

    /*
     * create square root matrix
     */
    A = MatrixXd(7, 7);
    A = P_arg.llt().matrixL();

    /*
     * create augmented sigma points
     */
    Xsig_aug.col(0) =  x_arg;

    for (i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_arg + sqrt(lambda_ + n_aug_) * A.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_arg - 
                                        sqrt(lambda_ + n_aug_) * A.col(i);
    }

    /********************************************************************
     * 1. Predict sigma points
     *
     ********************************************************************/

    for (i = 0; i < 2*n_aug_ + 1; i++) {

        /*
         * extract values for better readability
         */
        p_x = Xsig_aug(0,i);
        p_y = Xsig_aug(1,i);
        v = Xsig_aug(2,i);
        yaw = Xsig_aug(3,i);
        yawd = Xsig_aug(4,i);
        nu_a = Xsig_aug(5,i);
        nu_yawdd = Xsig_aug(6,i);

        /*
         * avoid division by zero
         */
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        } else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        v_p = v;
        yaw_p = yaw + yawd*delta_t;
        yawd_p = yawd;

        /*
         * add process noise longitudinal acceleration
         */
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        /*
         * add process noise yaw acceleration
         */
        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        /*
         * write predicted sigma point into right column
         */
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    /******************************************************************
     * 2. Predict state mean
     *
     ******************************************************************/
    
    /*
     * set weights
     */
    weights_(0) = lambda_/(lambda_ + n_aug_);

    for (i = 1; i < 2*n_aug_ + 1; i++) {
        weights_(i) = 0.5/(n_aug_ + lambda_);
    }

    x_.fill(0.0);
    for (i = 0; i < 2*n_aug_ + 1; i++) {
        x_ = x_ + weights_(i)*Xsig_pred_.col(i);
    }
    
    /**********************************************************************
     * 3. Predict state Covariance
     *
     * ********************************************************************/
    P_.fill(0.0);
    for (i = 0; i < 2*n_aug_ + 1; i++) {
    
        /*
         * state difference
         */
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        /*
         * angle normalization
         */
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ; 
    }

}

/*
 * UpdateLidar
 *
 * Updates the state and the state covariance matrix using a laser
 * measurement.
 * @in {MeasurementPackage} meas_package: has the current measurement
 *     for px and py from lidar
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  
    /*
     * Use lidar data to update the belief about the object's
     * position. Modify the state vector, x_, and covariance, P_.
     * Also calculate the lidar NIS.
     */

    /*
     * measurement dimension, px, py
     */
    int n_z = 2, i;
    double px, py;

    /*
     * Complete this function! Use radar data to update the belief about the
     * object's position. Modify the state vector, x_, and covariance, P_.
     */
     
    /*
     * 1. predict the sigmapoint in measurement space.
     */
    
    //create matrix for sigma points in measurement space
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd xyv = Xsig_pred_.col(i);
        px = xyv(0);
        py = xyv(1);
        Zsig.col(i) << px, py;
    }

    /*
     * 2. predict the mean in measurement space
     */
    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred += weights_(i)*Zsig.col(i);
    }

    /*
     * 3. predict the covariance matrix in measurement space
     */
    
    /*
     * measurement covariance matrix S
     */
    MatrixXd S = MatrixXd(n_z,n_z);

 
    S.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * diff * diff.transpose();
    }

    /*
     * Add Measurement Noise for Laser
     */
    S += R_laser;

    /*
     * UKF measurement update
     * 
     */
    px = meas_package.raw_measurements_[0];
    py = meas_package.raw_measurements_[1];
    VectorXd z = VectorXd(n_z);
    
    /*
     * real measurement data
     */
    z << px, py;

    /*
     * create matrix for cross correlation Tc
     */
    MatrixXd Tc = MatrixXd(n_x_, n_z); // 5x2

    /*
     * 1. calculate cross correlation matrix
     */
    Tc.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;       // 5 dimension 
        VectorXd z_diff = Zsig.col(i) - z_pred;         // 2 ""
        Tc += weights_(i)*x_diff*z_diff.transpose();    // 5x1 x 1x2 = 5x2
    }

    /*
     * 2. calculate Kalman gain K;
     */
    MatrixXd K = Tc*S.inverse();                        //5x2 x 2x2 = 5x2

    /*
     * 3. update state mean
     * x = 5x1, K = 5x2 , z -z_pred = 2x1
     */
    x_ = x_ + K*(z - z_pred);

    /*
     * 4.and covariance matrix
     * P = 5x5, K=5x2, S = 2x2 => 5x2*2x2 = (KS)5x2*2x5(K.transpose) = 5x5
     */
    P_ = P_ - K*S*K.transpose();

    /*
     * calculate the Lider NIS.
     */
    VectorXd m_diff = (z - z_pred);
    NIS_laser_ = m_diff.transpose()*S.inverse()*m_diff; // 1x2 x 2x2 x 2x1
                                                        // = 1x2 x 2x1 = 1x1    
}

/*
 * UpdateRadar
 *
 * Updates the state and the state covariance matrix using a radar
 * measurement.
 *
 * @in {MeasurementPackage} meas_package, contains 3 state
 *      measurement
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    
    /*
     * measurement dimension, radar can measure r, phi, and r_dot
     */
    int n_z = 3, i;
    double r, phi, r_dot, px, py, v, yaw;

    /*
     * Use radar data to update the belief about the object's position.
     * Modify the state vector, x_, and covariance, P_.
     */
     
    /*****************************************************************
     * 1. predict the sigmapoint in measurement space.
     *
     *****************************************************************/
    
    /*
     * Create matrix for sigma points in measurement space
     */
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd xyv = Xsig_pred_.col(i);
        px = xyv(0);
        py = xyv(1);
        v = xyv(2);
        yaw = xyv(3);

        r = sqrt(px*px + py*py);
        phi = atan2(py,px);
       
        /*
         * if range is close to zero then set the range rate to zero.
         * no change 
         */
        if (0 /* fabs(r) < 0.0001 */) {
            r_dot = 0;
        } else {
            r_dot = (px*cos(yaw)*v + py*sin(yaw)*v)/r;
        }
        Zsig.col(i) << r, phi, r_dot;
    }

    /******************************************************************
     * 2. predict the mean in measurement space
     *
     ******************************************************************/
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        z_pred += weights_(i)*Zsig.col(i);
    }

    /*******************************************************************
     * 3. predict the covariance matrix in measurement space
     * 
     ********************************************************************/
    
    /*
     * measurement covariance matrix S
     */
    MatrixXd S = MatrixXd(n_z,n_z);
 
    S.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd diff = Zsig.col(i) - z_pred;
            
            //angle normalization
        while (diff(1) > M_PI) diff(1) -= 2*M_PI;
        while (diff(1) < -M_PI) diff(1) += 2*M_PI;
        S = S + weights_(i) * diff * diff.transpose();
    }

    /*
     * Add Measurement Noise for Radar
     */
    S += R_radar;

    /*
     * UKF measurement update
     * 
     */
    r = meas_package.raw_measurements_[0];
    phi = meas_package.raw_measurements_[1];
    r_dot = meas_package.raw_measurements_[2];
    VectorXd z = VectorXd(n_z);

    /*
     * First time take the real measurement data
     */
    z << r, phi, r_dot;

    /*
     * create matrix for cross correlation Tc
     */
    MatrixXd Tc = MatrixXd(n_x_, n_z); // 5x3 matrix

    /*
     * 1. calculate cross correlation matrix
     */
    Tc.fill(0.0);
    for (i = 0; i < 2 * n_aug_ + 1; i++) {
        VectorXd x_diff = Xsig_pred_.col(i) - x_;       // 5 
        VectorXd z_diff = Zsig.col(i) - z_pred;         // 3
        Tc += weights_(i)*x_diff*z_diff.transpose();    // 5x1 x 1x3 = 5x3
    }

    /*
     * 2. calculate Kalman gain K;
     */
    MatrixXd K = Tc*S.inverse();                        //5x3 x 3x3 = 5x3

    /*
     * 3. update state mean
     * x = 5x1, K = 5x3 , z -z_pred = 3x1
     */
    x_ = x_ + K*(z - z_pred);

    /*
     * 4.and covariance matrix
     * P = 5x5, K=5x3, S = 3x3 => 5x3*3x3 = (KS)5x3*3x5(K.transpose) = 5x5
     */
    P_ = P_ - K*S*K.transpose();

    /*
     * Calculate the radar NIS.
     */
    VectorXd m_diff = (z - z_pred);
    NIS_radar_ = m_diff.transpose()*S.inverse()*m_diff; // 1x3 x 3x3 x 3x1
                                                        // = 1x3 x 3x1 = 1x1
    
}
