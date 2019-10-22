#include "kalman_filter.h"

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

void KalmanFilter::Predict()
{
  /**
   * predict the state
   */
  x_ = F_ * x_;                       // Predict state
  P_ = F_ * P_ * F_.transpose() + Q_; // Predict state covariance
}

void KalmanFilter::Update(const VectorXd &z)
{
  /**
   * update the state by using Kalman Filter equations
   */
  VectorXd err = z - H_ * x_;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // new state
  x_ += K * err;
  MatrixXd I_ = MatrixXd::Identity(4, 4);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * update the state by using Extended Kalman Filter equations
   */
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // Calculate nonlinear function h(x)
  VectorXd hx = VectorXd(z.size());
  float rho = sqrt(px*px + py*py);
  hx(0) = rho;  // rho
  hx(1) = atan2(py, px);  // phi
  hx(2) = (px*vx + py*vy)/rho;  // rho_dot
  // Calculate kalman error
  VectorXd err = z - hx; 

  // Update Kalman matrixes (note: H is the Jacobian matrix) 
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();

  // new state
  x_ += K * err;
  MatrixXd I_ = MatrixXd::Identity(4, 4);
  P_ = (I_ - K * H_) * P_;


}
