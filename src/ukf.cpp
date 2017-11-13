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
  P_ = MatrixXd::Identity(5, 5); 
  //std::cout << P_ << std::endl;
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.;
  //double std_a = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;
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
  n_sig_ = 2 * n_aug_;

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_sig_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, n_sig_ + 1);

  //define spreading parameter
  lambda_ = 3 - n_aug_;
  weights_ = VectorXd(n_sig_ + 1);
  //set weights
  weights_[0] = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_sig_ + 1; i++)
  {
	  double weight = 0.5 / (n_aug_ + lambda_);
	  weights_(i) = weight;
  }
  radar_n_z_ = 3;
  laser_n_z_ = 2;
  radar_NIS_ = 0;
  lidar_NIS_ = 0;

  radar_Zsig_ = MatrixXd(radar_n_z_, n_sig_ + 1);
  laser_Zsig_ = MatrixXd(laser_n_z_, n_sig_ + 1);
  radar_z_pred_= VectorXd(radar_n_z_);
  laser_z_pred_ = VectorXd(laser_n_z_);
  VectorXd laser_z_pred_;
  radar_S_ = MatrixXd(radar_n_z_, radar_n_z_);
  laser_S_ = MatrixXd(laser_n_z_, laser_n_z_);
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
			//vx = rho_dot * cos(phi);
			//vy = rho_dot * sin(phi);
			//v = sqrt(vx*vx + vy*vy);
			//yaw = atan2(vy, vx);
			x_ << px, py, 0, 0,0;
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
	if (abs(meas_package.raw_measurements_[0]) < 0.000001 || abs(meas_package.raw_measurements_[1]) < 0.000001) {
		return;
	}
	float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
	previous_timestamp_ = meas_package.timestamp_;

	/*****************************************************************************
	*  Prediction
	****************************************************************************/
	if (dt > 0.0001)
	{/*
		cout << "Prediction" << endl;*/
		Prediction(dt);
	}
	/*****************************************************************************
	*  Update
	****************************************************************************/
	/**
	TODO:
	* Use the sensor type to perform the update step.
	* Update the state and covariance matrices.
	*/
	//cout << "update" << endl;
	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		// Radar updates
		PredictRadarMeasurement();
		UpdateRadar(meas_package);
		//cout << "radar updated" << endl;
	}
	else {
		PredictLaserMeasurement();
		UpdateLidar(meas_package);/*
		cout << "lidar updated" << endl;*/
	}


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
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);
	PredictMeanAndCovariance();
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
	cout << "UpdateLidar" << endl;
	VectorXd z = VectorXd(laser_n_z_);
	z <<
		meas_package.raw_measurements_[0],   //px
		meas_package.raw_measurements_[1];   //py
	MatrixXd Tc = MatrixXd(n_x_, laser_n_z_);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

												//residual
		VectorXd z_diff = laser_Zsig_.col(i) - laser_z_pred_;
		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	MatrixXd K = Tc * laser_S_.inverse();
	VectorXd z_diff = z - laser_z_pred_;
	lidar_NIS_ = z_diff.transpose()*laser_S_.inverse()*z_diff;
	cout << "lidar_NIS_" << lidar_NIS_ <<endl;
	x_ = x_ + K * z_diff;
	P_ = P_ - K*laser_S_*K.transpose();
	
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
	VectorXd z = VectorXd(radar_n_z_);
	z <<
		meas_package.raw_measurements_[0],   //rho in m
		meas_package.raw_measurements_[1],   //phi in rad
		meas_package.raw_measurements_[2];   //rho_dot in m/s

	MatrixXd Tc = MatrixXd(n_x_, radar_n_z_);
	Tc.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		VectorXd z_diff = radar_Zsig_.col(i) - radar_z_pred_;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
	}
	MatrixXd K = Tc * radar_S_.inverse();
	VectorXd z_diff = z - radar_z_pred_;
	radar_NIS_ = z_diff.transpose()*radar_S_.inverse()*z_diff;
	cout << "radar_NIS_ " << radar_NIS_ <<endl;
	while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;
	//update 
	x_ = x_ + K * z_diff;
	P_ = P_ - K*radar_S_*K.transpose();
}

void UKF::AugmentedSigmaPoints() {


	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug_);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

	//create augmented mean state
	x_aug.head(5) = x_;
	//std::cout << "x = " << std::endl << x_ << std::endl;
	//std::cout << "x_aug = " << std::endl << x_aug << std::endl;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	//std::cout << "P_aug = " << std::endl << P_aug << std::endl;
	P_aug.topLeftCorner(5, 5) = P_;
	P_aug(5, 5) = std_a_*std_a_;
	P_aug(6, 6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//write augmented sigma points
	Xsig_aug_.col(0) = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug_.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug_.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
	}
}
void UKF::SigmaPointPrediction(double delta_t)
{
	//std::cout << "SigmaPointPrediction "<< std::endl;
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug_(0, i);
		double p_y = Xsig_aug_(1, i);
		double v = Xsig_aug_(2, i);
		double yaw = Xsig_aug_(3, i);
		double yawd = Xsig_aug_(4, i);
		double nu_a = Xsig_aug_(5, i);
		double nu_yawdd = Xsig_aug_(6, i);
		//predicted state values
		double px_p, py_p;
		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
		}
		else { //straight line
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}
		double v_p = v;
		double yaw_p = yaw + yawd*delta_t;
		double yawd_p = yawd;
		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;
		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;
		//write predicted sigma point into right column
		Xsig_pred_(0, i) = px_p;
		Xsig_pred_(1, i) = py_p;
		Xsig_pred_(2, i) = v_p;
		Xsig_pred_(3, i) = yaw_p;
		Xsig_pred_(4, i) = yawd_p;
	}
}
void UKF::PredictMeanAndCovariance(void)
{
	//predict state mean
	x_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		x_ = x_ + weights_(i)*Xsig_pred_.col(i);
	}
	//predict state covariance matrix
	P_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		VectorXd x_diff = Xsig_pred_.col(i) - x_;
		while (x_diff(3)> M_PI) x_diff(3) -= 2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3) += 2.*M_PI;

		P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
	}
}
void UKF::PredictRadarMeasurement() { 
	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);
		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;

		radar_Zsig_(0, i) = sqrt(p_x*p_x + p_y*p_y);                        //r
		radar_Zsig_(1, i) = atan2(p_y, p_x);                                 //phi
		radar_Zsig_(2, i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
	}

	//calculate mean predicted measurement
	radar_z_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		radar_z_pred_ = radar_z_pred_ + weights_(i) * radar_Zsig_.col(i);
	}

	//calculate measurement covariance matrix S
	radar_S_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
											   //residual
		VectorXd z_diff = radar_Zsig_.col(i) - radar_z_pred_;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1) -= 2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1) += 2.*M_PI;

		radar_S_ = radar_S_ + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(radar_n_z_, radar_n_z_);
	R << std_radr_*std_radr_, 0, 0,
		0, std_radphi_*std_radphi_, 0,
		0, 0, std_radrd_*std_radrd_;
	radar_S_ = radar_S_ + R;


	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	//std::cout << "radar_z_pred_: " << std::endl << radar_z_pred_ << std::endl;
	//std::cout << "radar_S_: " << std::endl << radar_S_ << std::endl;
}
void UKF::PredictLaserMeasurement() {
	//cout << "PredictLaserMeasurement" << endl;
	//transform sigma points into measurement space
	laser_Zsig_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++)
	{
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);

		laser_Zsig_(0, i) = p_x;                        //px
		laser_Zsig_(1, i) = p_y;                        //py
		
	}

	//calculate mean predicted measurement
	laser_z_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {
		laser_z_pred_ = laser_z_pred_ + weights_(i) * laser_Zsig_.col(i);
	}

	//calculate measurement covariance matrix S
	laser_S_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points residual
		VectorXd z_diff = laser_Zsig_.col(i) - laser_z_pred_;
		laser_S_ = laser_S_ + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd(laser_n_z_, laser_n_z_);
	R << std_laspx_*std_laspx_, 0,
		0, std_laspy_*std_laspy_;
	laser_S_ = laser_S_ + R;


	///*******************************************************************************
	//* Student part end
	//******************************************************************************/

	//////print result
	//std::cout << "laser_z_pred_: " << std::endl << laser_z_pred_ << std::endl;
	//std::cout << "laser_S_: " << std::endl << laser_S_ << std::endl;
}
