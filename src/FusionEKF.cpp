#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

enum E_SensorType {LASER, RADAR, BOTH};

// Set sensor to be used
E_SensorType sensor = BOTH;  

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  // Px and Py measerument
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // Hj_ is initialized when used
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  float dt = static_cast<float>(measurement_pack.timestamp_ - previous_timestamp_)/1000000.0; 
  previous_timestamp_ = measurement_pack.timestamp_;

  float noise_ax = 9;
  float noise_ay = 9;
  float dt2 = dt*dt;
  float dt4 = dt2*dt2;  

  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    //cout << "EKF: init" << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout << "EKF: init radar" << endl;
      // Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      float ro = measurement_pack.raw_measurements_(0);
      float theta = measurement_pack.raw_measurements_(1);
      // float theta_dot = measurement_pack.raw_measurements_(2);
      float Px = cos(theta)*ro;
      float Py = sin(theta)*ro;
      ekf_.x_ << Px, Py, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      cout << "EKF: init Laser" << endl;
      // Initialize state.
      float Px = measurement_pack.raw_measurements_(0);
      float Py = measurement_pack.raw_measurements_(1);
      ekf_.x_ << Px, Py, 0, 0;
    }
    //cout << "EKF: init Q" << endl;
    // Initialize P to a large value
    dt = 0.1; // Set to fixed value for first measurement to avoid very high values
    dt2 = dt*dt;
    dt4 = dt2*dt2;

    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt4/4.0 * noise_ax, 0, (dt2*dt/2.0)*noise_ax, 0,
            0, dt4/4.0 * noise_ay, 0, dt2*dt/2.0 * noise_ay,
            (dt2*dt/2.0)*noise_ax, 0, dt2*noise_ax, 0,
            0, (dt2*dt/2.0)*noise_ay, 0, dt2*noise_ay;
    //cout << "EKF: init P" << endl;
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ = 10.0*ekf_.Q_;

    // Initialize F
    ekf_.F_ = MatrixXd(4, 4);

    // done initializing, no need to predict or update
    is_initialized_ = true;
/*
    cout << "\nx_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    cout << "dt_ = " << dt << endl;
    cout << "EKF: init done!" << endl;
    */
    return;
  }
   
  /**
   * Prediction
   */

  /**
   *  Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  //cout << "\nEKF: Predict" << endl;
  // Update F
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1; 

  // Set the process covariance matrix Q

  //cout << "Update Q" << endl;
  ekf_.Q_ << dt4/4.0 * noise_ax, 0, (dt2*dt/2.0)*noise_ax, 0,
            0, dt4/4.0 * noise_ay, 0, dt2*dt/2.0 * noise_ay,
            (dt2*dt/2.0)*noise_ax, 0, dt2*noise_ax, 0,
            0, (dt2*dt/2.0)*noise_ay, 0, dt2*noise_ay;
  ekf_.Predict();
/*
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
  cout << "dt_ = " << dt << endl;
  */
  /**
   * Update
   */

  /**
   * 
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //cout << "\nEKF: Update radar" << endl;
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    if (Hj_(0, 0) != -1 && Hj_(1, 1) != -1)
    {
      ekf_.H_ = Hj_;
    } // Else use previous Hj value
      ekf_.R_ = R_radar_;
      if ((sensor == BOTH) or (sensor == RADAR))
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    
  } else {
    //cout << "\nEKF: Update laser" << endl;
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    if ((sensor == BOTH) or (sensor == LASER))
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "\nx_ = \n" << ekf_.x_ << endl;
  cout << "P_ = \n" << ekf_.P_ << endl;
}
