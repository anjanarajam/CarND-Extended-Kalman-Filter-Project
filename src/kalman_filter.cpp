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

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
    /*Projection 4D space(predicted value) to 2D space(measurement space)*/
    VectorXd z_pred = H_ * x_;

    /* Measurement update compares the measured value from the sensor to the predicted value */
    VectorXd y = z - z_pred;

    /*Calculate the estimated value of state and covariance */
    CalculateEstimatedValue(y);

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
    VectorXd z_pred(3);    

    /* Get the predicted values */
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    /*Maunually convert to polar cordinates*/
    double rho = sqrt(px * px + py * py);
    double phi = atan2(py, px);
    double rho_dot = (px * vx + py * vy) / rho;

    /* Fill the vector */
    z_pred << rho, phi, rho_dot;

    /* Measurement update compares the measured value from the 
    sensor to the predicted value */
    VectorXd y = z - z_pred;

    /*Calculate the estimated value of state and covariance */
    CalculateEstimatedValue(y);
}

void KalmanFilter::CalculateEstimatedValue(const VectorXd& y)
{
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    /*Estimated state */
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    /* Estimated covariance */
    P_ = (I - K * H_) * P_;
}