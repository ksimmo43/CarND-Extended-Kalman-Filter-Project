#include <iostream>
#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}
/**
void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_laser_ = H_in;
  R_laser_ = R_in;
  Hj_ = H_in;
  R_radar_ = R_in;
  Q_ = Q_in;
}
*/
void KalmanFilter::Predict() {
  x_ = F_ * x_;
  //std::cout << "P: " << P_ << std::endl;
  P_ = F_ * P_ * F_.transpose() + Q_;
  //std::cout << "P: " << P_ << std::endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  //VectorXd z_pred = H_laser_ * x_;
	VectorXd y = z - H_laser_ * x_;
	//MatrixXd Ht = H_.transpose();
	MatrixXd S = H_laser_ * P_ * H_laser_.transpose() + R_laser_;
	//MatrixXd Si = S.inverse();
	//MatrixXd PHt = P_ * H_laser_.transpose();
	MatrixXd K = P_ * H_laser_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  //recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  //pre-compute terms of z_pred based on cart to polar coordinate conversion
  float c1 = sqrt((px*px)+(py*py));

  if(fabs(c1) < 0.0001){
		std::cout << "Possible Division by Zero, denominator set = to 0.0001" << std::endl;
    c1 += 0.0001;
	}

  float c2 = atan2(py, px);
  float c3 = ((px*vx)+(py*vy))/c1;

  VectorXd z_pred(3);
  z_pred <<  c1,
        c2,
        c3;

	VectorXd y = z - z_pred;

  //Make sure phi is  -pi < phi < pi
  y(1) = tools.normalizeAngle(y(1));
	//MatrixXd Ht = H_.transpose();
	MatrixXd S = Hj_ * P_ * Hj_.transpose() + R_radar_;
	//MatrixXd Si = S.inverse();
	//MatrixXd PHt = P_ * H_laser_.transpose();
	MatrixXd K = P_ * Hj_.transpose() * S.inverse();

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * Hj_) * P_;
}
