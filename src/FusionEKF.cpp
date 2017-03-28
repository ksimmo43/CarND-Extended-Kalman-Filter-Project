#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  ekf_.x_ = VectorXd(4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.R_laser_ = MatrixXd(2, 2);
  ekf_.R_radar_ = MatrixXd(3, 3);
  ekf_.H_laser_ = MatrixXd(2, 4);
  ekf_.Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  ekf_.R_laser_ << 0.0225, 0,
        0, 0.0225;
  //measurement covariance matrix - radar
  ekf_.R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  //measurement matrix - laser
  ekf_.H_laser_ << 1, 0, 0, 0,
        0, 1, 0, 0;
  //state covariance matrix P
  ekf_.P_ << 1, 0, 0, 0,
  		  0, 1, 0, 0,
  		  0, 0, 1000, 0,
  		  0, 0, 0, 1000;
  //the initial transition matrix F_
  ekf_.F_ << 1, 0, 1, 0,
  		  0, 1, 0, 1,
  		  0, 0, 1, 0,
  		  0, 0, 0, 1;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    // first measurement
    cout << "EKF: " << endl;
    cout << "P: " << ekf_.P_ << endl;
    cout << "F: " << ekf_.F_ << endl;

    cout << "Rl: " << ekf_.R_laser_ << endl;
    cout << "Rr: " << ekf_.R_radar_ << endl;
    cout << "H: " << ekf_.H_laser_ << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho_i = measurement_pack.raw_measurements_(0);
      float phi_i = measurement_pack.raw_measurements_(1);

      // Make sure phi is  -pi < phi < pi
      phi_i = tools.normalizeAngle(phi_i);

      ekf_.x_ << rho_i * cos(phi_i), rho_i * sin(phi_i), 0, 0;
      previous_timestamp_ = measurement_pack.timestamp_;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

      previous_timestamp_ = measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;

    //define pi
    const double pi = 2*acos(0.0);

    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

   //compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;


  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  //std::cout << "F_= " << kf_.F_ << std::endl;
  //set acceleration noises
  float noise_ax = 9;
  float noise_ay = 9;
  //set the process covariance matrix Q
  //kf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
        0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
        0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // calculate Jacobian at current predicted state
    ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);

    //measurement update
  	ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    //measurement update
  	ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
