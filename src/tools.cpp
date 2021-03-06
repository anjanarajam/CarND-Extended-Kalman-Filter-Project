#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

/**
 * Constructor.
 */
Tools::Tools() {}

/**
 * Destructor.
 */
Tools::~Tools() {}

/**
 * Calculate the RMSE here.
 */
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	/** Check the validity of the following inputs:
	 * the estimation vector size should not be zero
	 * the estimation vector size should equal ground truth vector size
	 */
	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		std:: cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	/* Accumulate squared residuals */
	for (unsigned int i = 0; i < estimations.size(); ++i) {

		VectorXd residual = estimations[i] - ground_truth[i];

		/* Coefficient-wise multiplication */
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	/* Calculate the mean */
	rmse = rmse / estimations.size();

	/* Calculate the squared root */
	rmse = rmse.array().sqrt();

	/* Return the result */
	return rmse;
}

/**
 * Calculate the Jacobian matrix for linear approximation of 
 * h(x') in case of Radar.
 */
MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

	/* Initialze the measurement matrix for Radar */
	MatrixXd Hj(3, 4);

	// Get the state parameters
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);

	/* Calculate denominators */
	double denom_1 = px * px + py * py;
	double denom_2 = sqrt(denom_1);
	double denom_3 = (denom_1 * denom_2);

	/* Check division by zero error */
	if (fabs(denom_1) < 0.0001) {
		std:: cout << "ERROR:Division by Zero" << std::endl;
		return Hj;
	}

	/* Calculate Jacobian matrix */
	Hj << (px / denom_2), (py / denom_2), 0, 0,
		-(py / denom_1), (px / denom_1), 0, 0,
		py* (vx * py - vy * px) / denom_3, px* (px * vy - py * vx) / denom_3, px / denom_2, py / denom_2;

	return Hj;
}
