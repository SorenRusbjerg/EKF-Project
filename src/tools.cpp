#include "tools.h"
#include <iostream>
#include <math.h>       /* pow */

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size()==0 or estimations.size() != ground_truth.size()){
    std::cout << "Size error" << std::endl;
    return rmse;  
  }
  VectorXd res(4);
  res << 0,0,0,0;
  for (unsigned int i=0; i < estimations.size(); ++i) {
    // ... your code here
    VectorXd sq = estimations[i] - ground_truth[i];
    sq = sq.array() * sq.array();
    res += sq;
    
  }

  // calculate the mean
  res = res/estimations.size();

  // calculate the squared root
  rmse = res.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  // check division by zero
  if (px==0.0 and py==0.0)
  {
      std::cout << "division by zero" << std::endl;
      return Hj;
  }
  
  // compute the Jacobian matrix
  Hj << px/sqrt(px*px + py*py), py/sqrt(px*px + py*py), 0.0, 0.0,
        -py/(px*px + py*py), px/(px*px + py*py), 0.0, 0.0,
        py*(vx*py-vy*px)/pow(px*px + py*py, 1.5), px*(vy*px-vx*py)/pow(px*px + py*py, 1.5), px/sqrt(px*px + py*py), py/sqrt(px*px + py*py);

  return Hj;
}
