#include "ukf.h"
#include "Eigen/Dense"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_.fill(0);
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 8.7;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.695;

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
   n_x_ = 5;

   n_aug_ = 7;

   n_sig_points  = 2 * n_aug_ + 1;

   weights_ = VectorXd(n_sig_points);

   weights_.fill(0);

   Xsig_pred_ = MatrixXd(n_x_, n_sig_points);

   Xsig_pred_.fill(0);

   lambda_ = 3 - n_aug_;

   is_initialized_ = false;

   InitializeWeights();
}

UKF::~UKF() {}

void UKF::InitializeWeights() {
    weights_.fill(0);
    weights_[0] = lambda_/(lambda_+n_aug_);
    for(int i {1}; i < 2*n_aug_+1; ++i){
        weights_[i] = 0.5/(lambda_+n_aug_);
    }
}

void UKF::ProcessMeasurement(MeasurementPackage &meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   if(!is_initialized_){
       if((meas_package.sensor_->sensor_type_ == Sensor::LASER) && use_laser_){
           x_[0] = meas_package.raw_measurements_[0];
           x_[1] = meas_package.raw_measurements_[1];
           x_[2]  = 0;
           x_[3]  = 0;
           x_[4]  = 0;

           P_ << pow(LidarSensor::std_laspx_,2), 0, 0, 0, 0,
                   0, pow(LidarSensor::std_laspy_,2), 0, 0, 0,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0;

           is_initialized_ = true;
       }else if((meas_package.sensor_->sensor_type_ == Sensor::RADAR) && use_radar_){
           x_[0] = meas_package.raw_measurements_[0] * cos(meas_package.raw_measurements_[1]);
           x_[1] = meas_package.raw_measurements_[0] * sin(meas_package.raw_measurements_[1]);
           x_[2]  = 0;
           x_[3]  = 0;
           x_[4]  = 0;

           P_ << RadarSensor::std_radr_*cos(RadarSensor::std_radphi_), 0, 0, 0, 0,
                   0, RadarSensor::std_radr_*sin(RadarSensor::std_radphi_), 0, 0, 0,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0,
                   0, 0, 0, 0, 0;

           is_initialized_ = true;
       }
       time_us_ = meas_package.timestamp_;
       return;
   }

    const double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);

    if((meas_package.sensor_->sensor_type_ == meas_package.sensor_->RADAR) && use_radar_){

        Update(meas_package, nis_radar);
    }else if((meas_package.sensor_->sensor_type_ == meas_package.sensor_->LASER) && use_laser_){
        Update(meas_package, nis_laser);
    }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
    Eigen::VectorXd x_aug_{Eigen::VectorXd::Zero(n_aug_)};
    Eigen::MatrixXd P_aug_ {Eigen::MatrixXd::Zero(n_aug_, n_aug_)};
    Eigen::MatrixXd Q_{MatrixXd::Zero(2,2)};
    MatrixXd Xsig_aug_ {MatrixXd(n_aug_, 2*n_aug_+1)};

    x_aug_.segment(0,5) = x_;
    x_aug_[5] = 0;
    x_aug_[6] = 0;

    Q_ << pow(std_a_,2) , 0, 0, pow(std_yawdd_,2);
    P_aug_.block<5,5>(0,0) = P_;
    P_aug_.block<2,2>(5,5) = Q_;

    MatrixXd sqrt_P_  = P_aug_.llt().matrixL();

    Xsig_aug_.col(0) = x_aug_;
    for(int i {0}; i < n_aug_; ++i){
        Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_+n_aug_) * sqrt_P_.col(i);
        Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * sqrt_P_.col(i);
    }

    //Prediction function
    for(int i {0}; i < 2 * n_aug_ +1; ++i){
        // extract values for better readability
        const double p_x = Xsig_aug_(0,i);
        const double p_y = Xsig_aug_(1,i);
        const double v = Xsig_aug_(2,i);
        const double yaw = Xsig_aug_(3,i);
        const double yawd = Xsig_aug_(4,i);
        const double nu_a = Xsig_aug_(5,i);
        const double nu_yawdd = Xsig_aug_(6,i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( - cos(yaw+yawd*delta_t) + cos(yaw));
        } else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        // add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        // write predicted sigma point into right column
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
    }

    x_.fill(0.0);
    x_ = Xsig_pred_ *  weights_;

    P_.fill(0.0);
    for(int i {0}; i < 2*n_aug_+1; ++i){
        VectorXd x_diff {Xsig_pred_.col(i) - x_};

        verify_state_parameter(x_diff);
        P_ = P_ + weights_[i] * (x_diff*x_diff.transpose());
    }
}

void UKF::Update(MeasurementPackage &meas_package, std::vector<double> &nis) {
    /**
     * TODO: Complete this function! Use radar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the radar NIS, if desired.
     */
    int n_meas_ = meas_package.sensor_->n_meas_;

    MatrixXd Z_sig_pkt {MatrixXd::Zero(n_meas_, n_sig_points)};
    MatrixXd S_ {MatrixXd::Zero(n_meas_, n_meas_)};
    MatrixXd z_diff {MatrixXd::Zero(n_meas_, n_sig_points)};
    VectorXd z_ {VectorXd::Zero(n_meas_)};
    MatrixXd T_ {MatrixXd::Zero(n_x_, n_meas_)};

    meas_package.sensor_->measurement_function(Xsig_pred_, Z_sig_pkt, n_sig_points);

    for (int i{0}; i < n_sig_points; ++i) {
        z_ = z_ + weights_[i] * Z_sig_pkt.col(i);
        z_diff.col(i) = Z_sig_pkt.col(i) - z_;

        VectorXd tmp_vector {z_diff.col(i)};
        meas_package.sensor_->verify_parameter(tmp_vector);

        S_ = S_ + weights_[i] * z_diff.col(i) * z_diff.col(i).transpose();
    }

    S_ = S_ + meas_package.sensor_->r_matrix();;

    VectorXd z_diff_meas {meas_package.raw_measurements_ - z_};

    meas_package.sensor_->verify_parameter(z_diff_meas);

    for (int i{0}; i < n_sig_points; ++i) {
        VectorXd x_diff {Xsig_pred_.col(i) - x_};
        verify_state_parameter(x_diff);

        T_ = T_ +  weights_[i] * x_diff * z_diff.col(i).transpose();
    }

    MatrixXd kalman_gain{T_ * S_.inverse()};

    x_ = x_ + kalman_gain * z_diff_meas;
    P_ = P_ - kalman_gain * S_ * kalman_gain.transpose();

    nis.push_back(z_diff_meas.transpose()*S_.inverse()*z_diff_meas);
};

void UKF::verify_state_parameter(Eigen::VectorXd &x_diff){
    if (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    if (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
}