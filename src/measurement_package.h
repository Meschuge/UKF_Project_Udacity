#ifndef MEASUREMENT_PACKAGE_H_
#define MEASUREMENT_PACKAGE_H_

#include "Eigen/Dense"
class Sensor {
public:
    enum SensorType{
        LASER,
        RADAR
    } sensor_type_;

    int n_meas_ = -1;

    Sensor(int n_): n_meas_(n_) {}
    virtual void measurement_function(const Eigen::MatrixXd &x_predicted_sigma, Eigen::MatrixXd &z_sigma_point, const int &n_) = 0;
    virtual void verify_parameter(Eigen::VectorXd &z_) = 0;
    virtual Eigen::MatrixXd r_matrix() = 0;
};

class MeasurementPackage {
public:
    MeasurementPackage() {};
    ~MeasurementPackage() { delete sensor_;};
    long timestamp_;
    Eigen::VectorXd raw_measurements_;

    Sensor *sensor_ = nullptr;
};

class RadarSensor: public Sensor{
public:
    constexpr static double std_radr_ = 0.3;
    constexpr static double std_radphi_ = 0.03;
    constexpr static double std_radrd_ = 0.3;

    RadarSensor():Sensor(3){
        sensor_type_ = Sensor::RADAR;
    }
    void measurement_function(const Eigen::MatrixXd &x_predicted_sigma, Eigen::MatrixXd &z_sigma_point, const int &n_){
          for (int i{0}; i < n_; ++i) {
              const double p_x = x_predicted_sigma(0,i);
              const double p_y = x_predicted_sigma(1,i);
              const double v = x_predicted_sigma(2,i);
              const double yaw = x_predicted_sigma(3,i);

              const double v1 = cos(yaw)*v;
              const double v2 = sin(yaw)*v;

              z_sigma_point(0,i) = sqrt(p_x*p_x + p_y*p_y);
              z_sigma_point(1,i) = atan2(p_y, p_x);
              z_sigma_point(2,i) = (p_x*v1 + p_y*v2)/sqrt(p_x*p_x + p_y*p_y);
          }
    }
    void verify_parameter(Eigen::VectorXd &z_){
        if(z_(1) > M_PI) z_(1) -= 2.*M_PI;
        if(z_(1) < -M_PI) z_(1) += 2.*M_PI;
    }
    Eigen::MatrixXd r_matrix(){
        Eigen::MatrixXd R_ = Eigen::MatrixXd::Zero(3, 3);
        R_ << pow(std_radr_, 2), 0, 0, 0, pow(std_radphi_, 2), 0, 0, 0, pow(std_radrd_, 2);

        return R_;
    }
};

class LidarSensor: public Sensor{
public:
    constexpr static double std_laspx_ = 0.35;
    constexpr static double std_laspy_ = 0.15;

    LidarSensor():Sensor(2){
        sensor_type_ = Sensor::LASER;
    }
    void measurement_function(const Eigen::MatrixXd &x_predicted_sigma, Eigen::MatrixXd &z_sigma_point, const int & n_){
        for (int i {0}; i < n_ ; ++i) {
            z_sigma_point(0, i) = x_predicted_sigma(0,i);
            z_sigma_point(1, i) = x_predicted_sigma(1,i);
        }
    }
    void verify_parameter(Eigen::VectorXd &z_){}
    Eigen::MatrixXd r_matrix(){
        Eigen::MatrixXd R_ = Eigen::MatrixXd::Zero(2, 2);
        R_ << pow(std_laspx_, 2), 0, 0, pow(std_laspy_, 2);
        return R_;
    }

};

#endif /* MEASUREMENT_PACKAGE_H_ */