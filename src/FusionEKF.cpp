#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  /* initializing matrices */
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  /* measurement covariance matrix for laser */
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  /* measurement covariance matrix for radar */
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  /* Set the state transition matrix */ 
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
      0, 1, 0, 1,
      0, 0, 1, 0,
      0, 0, 0, 1;

  /* Set state covariance matrix P */
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

/* Set the measurement matrix for laser*/
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
    /* If not yet initialized */
    if (!is_initialized_) {
        /**
         * TODO: Initialize the state ekf_.x_ with the first measurement.
         * TODO: Create the covariance matrix.
         * You'll need to convert radar from polar to cartesian coordinates.
         */
        /* Initialize position and velocity */
        double px{};
        double py{};
        double vx{};
        double vy{};
    
        // first measurement 
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
          // TODO: Convert radar from polar to cartesian coordinates 
          //         and initialize state.
            /* Get range - radial distance from the origin */
            double rho = measurement_pack.raw_measurements_[0];
            /* Get bearing - angle between rho and state*/
            double phi = measurement_pack.raw_measurements_[1];
            /* Get radial velocity -  rate of change of rho */
            double rho_dot = measurement_pack.raw_measurements_[0];

            /* x- position */
            px = rho * cos(phi);
            /* y - position */
            py = rho * sin(phi);
            /* x-velocity */
            vx = rho_dot * cos(phi);
            /* y - velocity */
            vy = rho_dot * sin(phi);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
          // TODO: Initialize state.
            /* x- position */
            px = measurement_pack.raw_measurements_[0];
            /* y - position */
            py = measurement_pack.raw_measurements_[1];
            /* x-velocity */
            vx = 0;
            /* y - velocity */
            vy = 0;
        }

        /* Get the initial state of the pedestrian*/
        ekf_.x_ << px, py, vx, vy;

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
  }

    /**
    * Prediction
    */

    /**
    * TODO: Update the state transition matrix F according to the new elapsed time.
    * Time is measured in seconds.
    * TODO: Update the process noise covariance matrix.
    * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
    */
    /* Initialize acceleration noise components */
    double noise_ax = 9;
    double noise_ay = 9;

    /* compute the time elapsed between the current and previous measurements
    dt - expressed in seconds */
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = measurement_pack.timestamp_;

    /* Calculate multiples of elapsed time */
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    /* Modify the F matrix so that the time is integrated */
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    // set the process covariance matrix Q
    ekf_.Q_ = MatrixXd(4, 4);
    ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
        0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
        dt_3 / 2 * noise_ax, 0, dt_2* noise_ax, 0,
        0, dt_3 / 2 * noise_ay, 0, dt_2* noise_ay;

    ekf_.Predict();

    /**
    * Update
    */
    /**
    * TODO:
    * - Use the sensor type to perform the update step.
    * - Update the state and covariance matrices.
    */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
        /* raw_measurements_ contains rho, phi and rho_dot for radar */
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;

        ekf_.UpdateEKF(measurement_pack.raw_measurements_);

    } else {
    // TODO: Laser updates
        // measurement update
        Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.H_ = Hj_;
        ekf_.R_ = R_radar_;
        /* raw_measurements_ contains px and py for lidar */
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    /* Print the new estimates */
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
