#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define DEBUG (0)
#define THRESHOLD (0.00001)
#define EPSILON (0.0001)

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // protected:
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2; // CTRV model; n_x_ + 2 (linear acceleration and radial acceleration noise)

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* the current NIS for radar
  NIS_radar_ = 0; // 7.815 | 3-degrees-of-freedom. chi-squared distribution table

  ///* the current NIS for laser
  NIS_laser_ = 0; // 5.991 | 2-degrees-of-freedom. chi-squared distribution table

  ///* Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_/(lambda_+n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {
    weights_(i) = 0.5/(n_aug_+lambda_);
  }

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_ << -1, -1, -1, -1, -1;

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_ <<
      100, 0, 0, 0, 0,
      0, 100, 0, 0, 0,
      0, 0, 100, 0, 0,
      0, 0, 0, 100, 0,
      0, 0, 0, 0, 100;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.8;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // Initialize predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  /**
  TODO: Done

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off... (done)
  */


}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {
  /**
  TODO:
  Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = (float) measurement_pack.raw_measurements_[0];
      float phi = (float) measurement_pack.raw_measurements_[1];
      float rodot = (float) measurement_pack.raw_measurements_[2];
      // TODO(Manav): Optional. Try to incorporate rodot for this init
      x_ << ro * cos(phi), ro * sin(phi), 0, 0, 0;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0, 0;
    }

    time_us_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  double dt = (measurement_pack.timestamp_ - time_us_) / 1000000.0;
  time_us_ = measurement_pack.timestamp_;

  // Predict measurements considering history and delta t
  UKF::Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Update state and covariance matrix considering current measurement
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    if(!use_radar_) return;
    UpdateRadar(measurement_pack);
  }
  else {
    if(!use_laser_) return;
    UpdateLidar(measurement_pack);
  }

}

static inline void normalize_angle(double &angle)
{
  while (angle> M_PI) angle-=2.*M_PI;
  while (angle<-M_PI) angle+=2.*M_PI;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /*
   * TODO:
   * Estimate the object's location. Modify state vector, x_. Predict sigma points, the state, and
   * state covariance matrix.
   */

  /*  Steps involved: */
  /* 1. Augment state and covariance matrix to consider process noise */
  /* 2. Compute augmented sigma points */
  /* 3. Predict sigma points */
  /* 4. Predict state mean */
  /* 5. Predict state covariance matrix */
  /**/

  /* 1. Augment state and covariance matrix to consider process noise */
  // Ref: Lesson 7.18: Augmentation Assignment
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);

  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  /* 2. Compute (augmented) sigma points */
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Ref: Lesson 7.15: Generating Sigma Points Assignment
  // set first column of sigma point matrix
  Xsig_aug.col(0) = x_aug;
  //set remaining sigma points
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  /* 3. Predict sigma points */
  // Ref: Lesson 7.21: Generating Sigma Points Assignment
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > THRESHOLD) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  /* 4. Predict state mean */
  // Ref: Lesson 7.24: Predicted Mean and Covariance Assignment
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  /* 5. Predict state covariance matrix */
  // Ref: Lesson 7.24: Predicted Mean and Covariance Assignment
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    normalize_angle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /**  Steps involved: */
  /* 1. Transform sigma points from process space to measurement space */
  /* 2. Compute mean prediction measurement */
  /* 3. Compute measurement covariance matrix */
  /* 4. Compute cross correlation matrix */
  /* 5. Update state mean and covariance matrix */
  /* 6. Compute NIS value */

  int n_z = 2;
  VectorXd z = VectorXd(2);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);

  /* 1. Transform sigma points from process space to measurement space */
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;
    Zsig(1,i) = p_y;
  }

  /* 2. Compute mean prediction measurement */
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  /* 3. Compute measurement covariance matrix */
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd Zsig_diff = MatrixXd(n_z, 2*n_aug_+1);

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    Zsig_diff.col(i) = Zsig.col(i) - z_pred; //residual
    normalize_angle(Zsig_diff.col(i)(1)); //angle normalization

    S = S + weights_(i) * Zsig_diff.col(i) * Zsig_diff.col(i).transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  S = S + R;

  /* 4. Compute cross correlation matrix */
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    VectorXd x_diff = Xsig_pred_.col(i) - x_; // state difference
    normalize_angle(x_diff(3)); //angle normalization

    Tc = Tc + weights_(i) * x_diff * Zsig_diff.col(i).transpose();
  }

  MatrixXd K = Tc * S.inverse(); //Kalman gain K;

  VectorXd z_diff = z - z_pred; //residual

  /* 5. Update state mean and covariance matrix */
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /* 6. Compute NIS value */
  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  if(DEBUG) cout << "NISL: " << NIS_laser_ << endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  /**  Steps involved: */
  /* 1. Transform sigma points from process space to measurement space */
  /* 2. Compute mean prediction measurement */
  /* 3. Compute measurement covariance matrix */
  /* 4. Compute cross correlation matrix */
  /* 5. Update state mean and covariance matrix */
  /* 6. Compute NIS value */

  int n_z = 3;
  VectorXd z = VectorXd(3);
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  z(0) = meas_package.raw_measurements_(0);
  z(1) = meas_package.raw_measurements_(1);
  z(2) = meas_package.raw_measurements_(2);

  /* 1. Transform sigma points from process space to measurement space */
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  /* 2. Compute mean prediction measurement */
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  /* 3. Compute measurement covariance matrix */
  MatrixXd S = MatrixXd(n_z,n_z);
  MatrixXd Zsig_diff = MatrixXd(n_z, 2*n_aug_+1);

  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    Zsig_diff.col(i) = Zsig.col(i) - z_pred; //residual

    normalize_angle(Zsig_diff.col(i)(1));  //angle normalization

    S = S + weights_(i) * Zsig_diff.col(i) * Zsig_diff.col(i).transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0, std_radrd_*std_radrd_;
  S = S + R;

  /* 4. Compute cross correlation matrix */
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    VectorXd x_diff = Xsig_pred_.col(i) - x_; // state difference

    normalize_angle(x_diff(3)); //angle normalization

    Tc = Tc + weights_(i) * x_diff * Zsig_diff.col(i).transpose();
  }

  MatrixXd K = Tc * S.inverse(); //Kalman gain K;

  VectorXd z_diff = z - z_pred; //residual

  normalize_angle(z_diff(1)); //angle normalization

  /* 5. Update state mean and covariance matrix */
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  /* 6. Compute NIS value */
  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  if(DEBUG) cout << "NISR: " << NIS_radar_ - 7.8 << endl;
}

#pragma clang diagnostic pop