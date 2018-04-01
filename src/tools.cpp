#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
    // Calculate the RMSE here.
    // The code written in this file follows the lab exercises written in Udacity lectures
    // and additional remarks from their presented material.

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    // check the validity of the following inputs:
    //  * the estimation vector size should not be zero
    //  * the estimation vector size should equal ground truth vector size
    if(estimations.size() != ground_truth.size() || estimations.size() == 0){
        cout << "Invalid estimation or ground_truth data" << endl;
        return rmse;
    }

    //accumulate squared residuals
    for(int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];

        //coefficient-wise multiplication
        residual = residual.array()*residual.array();
        rmse += residual;
    }

    //calculate the mean, then
    //calculate the squared root
    rmse = (rmse/estimations.size()).array().sqrt();

    //return the result
    return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    //Calculate a Jacobian here.

    //recover state parameters
    double px = x_state(0);
    double py = x_state(1);
    double vx = x_state(2);
    double vy = x_state(3);

    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px * px + py * py;
    double c2 = sqrt(c1);
    double c3 = (c1 * c2);

    Eigen::MatrixXd jacobian(3, 4);
    jacobian << (px / c2), (py / c2), 0, 0,
            -(py / c1), (px / c1), 0, 0,
            py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

    return jacobian;
}
