#ifndef UKF_H
#define UKF_H

#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage &meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a radar or lidar measurement
   * @param meas_package The measurement at k+1
   */
  void Update(MeasurementPackage &meas_package, std::vector<double> &nis);

  /**
   * Initialize weights for prediction and update steps.
   */
   void InitializeWeights();

   /**
    * Update the state x_ and the covariance matrice P_ based on the sensor specific calculations
    */
    void UpdateMeasurements(Eigen::MatrixXd &S_, Eigen::MatrixXd &z_diff, Eigen::VectorXd &z_diff_meas, int &n_,  std::vector<double> &nis);

    /**
     * Verify if the angle of the state parameter are inbound.
     */
     void verify_state_parameter(Eigen::VectorXd &x_diff);
  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  //Number of Sigmapoints
  int n_sig_points;

  // Sigma point spreading parameter
  double lambda_;

  std::vector<double> nis_radar;

  std::vector<double> nis_laser;

};

#endif  // UKF_H