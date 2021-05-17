#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {  
  P_ = F_ * P_ * F_.transpose() + Q_;
	x_ = F_ * x_;
}

void KalmanFilter::Update(const VectorXd &z) {
    
    VectorXd y = z - H_ * x_;
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();
    
    long size = x_.size();
    MatrixXd I = MatrixXd::Identity(size, size);

    // new state
    P_ = ( I - ( K * H_ ) ) * P_;
    x_ = x_ + ( K * y );
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float c = sqrt( (px * px) + ( py * py ) );

  VectorXd h1 = VectorXd(3);
  h1 << c, atan2(py, px), (px*vx + py*vy) / c;

  VectorXd y = z - h1;
  MatrixXd Ht = H_.transpose();

  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  long size = x_.size();
  MatrixXd I = MatrixXd::Identity(size, size);
 
  // calculate new state
  P_ = ( I - ( K * H_ ) ) * P_;
  x_ = x_ + ( K * y );
}
