#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
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
  P_ = F_ * P_ * F_.transpose();
	x_ = (F_ * x_) + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    
    VectorXd y = z - H_ * x_;
    VectorXd S = H_ * P_ * H_.transpose() + R_;
    VectorXd K = P_ * H_.transpose() * S.inverse();
    
    long x_size = x_.size();
    VectorXd I = MatrixXd::Identity(2, 2);

    // new state
    P_ = ( I - ( K * H_ ) ) * P_;
    x_ = x_ + ( K * y );
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  float c = sqrt( (px * px) + ( py * py ) );

  VectorXd h1 = VectorXd(3);
  h1 << c, atan2(py, px), (px*vx + py*vy) / c, atan2(py, px);

  VectorXd y = z - h1;
  VectorXd S = H_ * P_ * H_.transpose() + R_;
  VectorXd K = P_ * H_.transpose() * S.inverse();
  
  long x_size = x_.size();
  VectorXd I = MatrixXd::Identity(2, 2);

  // new state
  P_ = ( I - ( K * H_ ) ) * P_;
  x_ = x_ + ( K * y );
}
