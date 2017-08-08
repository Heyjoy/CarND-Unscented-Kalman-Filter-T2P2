#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(estimations[0].size());
	VectorXd sum(estimations[0].size());

	sum.fill(0);
	for (int i = 0; i<estimations.size(); i++)
	{
		VectorXd error = estimations[i] - ground_truth[i];
		VectorXd error2 = error.array()*error.array();
		sum += error2;
	}
	rmse = (sum / estimations.size()).array().sqrt();

	//cout <<"size is :" <<estimations.size()<< ground_truth.size() << endl;
	return rmse;
}