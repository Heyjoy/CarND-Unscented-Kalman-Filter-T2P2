#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.fill(0.0);
  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0,
	  0, 0, 1000, 0, 0,
	  0, 0, 0, 1000, 0,
	  0, 0, 0, 0, 1000;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;
  //double std_a = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  //double std_yawdd = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //define spreading parameter
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(2 * n_aug_ + 1);
  //set weights
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++)
  {
	  double weight = 0.5 / (n_aug_ + lambda_);
	  weights_(i) = weight;
  }
  //for the first measurment
  is_initialized_ = false;
  // previous timestamp
  previous_timestamp_=0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
	if (use_laser_ == false && meas_package.sensor_type_ == MeasurementPackage::LASER) return;
	if (use_radar_ == false && meas_package.sensor_type_ == MeasurementPackage::RADAR) return;
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
	if (!is_initialized_) {
		// first measurement
		cout << "UKF:initializ " << endl;
		float  px, py, v,vx, vy, yaw;

		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/
			float rho = meas_package.raw_measurements_[0];
			float phi = meas_package.raw_measurements_[1];
			float rho_dot = meas_package.raw_measurements_[2];
			px = rho * cos(phi);
			py = rho * sin(phi);
			vx = rho_dot * cos(phi);
			vy = rho_dot * sin(phi);
			v = sqrt(vx*vx + vy*vy);
			yaw = atan2(vy, vx);
			x_ << px, py, v, yaw,0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
			/**
			Initialize state.
			*/
			px = meas_package.raw_measurements_(0);
			py = meas_package.raw_measurements_(1);
			x_ << px, py, 0, 0, 0;
		}

		// done initializing, no need to predict or update
		is_initialized_ = true;
		previous_timestamp_ = meas_package.timestamp_;
		return;
	}
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = meas_package.timestamp_;


}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

	//set state dimension
	int n_x = 5;

	//define spreading parameter
	double lambda = 3 - n_x;

	//set example state
	VectorXd x = VectorXd(n_x);
	x << 5.7441,
		1.3800,
		2.2049,
		0.5015,
		0.3528;

	//set example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
		-0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
		0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
		-0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
		-0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

	//create sigma point matrix
	MatrixXd Xsig = MatrixXd(n_x, 2 * n_x + 1);

	//calculate square root of P
	MatrixXd A = P.llt().matrixL();

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//set first column of sigma point matrix
	Xsig.col(0) = x;

	//set remaining sigma points
	for (int i = 0; i < n_x; i++)
	{
		Xsig.col(i + 1) = x + sqrt(lambda + n_x) * A.col(i);
		Xsig.col(i + 1 + n_x) = x - sqrt(lambda + n_x) * A.col(i);
	}

	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	//std::cout << "Xsig = " << std::endl << Xsig << std::endl;

	//write result
	*Xsig_out = Xsig;

	/* expected result:
	Xsig =
	5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441
	1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38
	2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049
	0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015
	0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879
	*/

}