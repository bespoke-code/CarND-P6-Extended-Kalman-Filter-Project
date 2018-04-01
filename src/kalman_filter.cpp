#include "kalman_filter.h"
#include <iostream>
#include <cmath>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {

    x_ = VectorXd(4);
    x_ << 0.0, 0.0, 0.0, 0.0;
    //state covariance matrix P
    P_ = MatrixXd(4, 4);
    P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;


    //measurement covariance
    //R_ = MatrixXd(2, 2);
    //R_ << 0.0225, 0,
    //        0, 0.0225;

    //measurement matrix
    H_ = MatrixXd(2, 4);
    H_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    //the initial transition matrix F_
    F_ = MatrixXd(4, 4);
    F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    Q_ = MatrixXd(4, 4);

    process_noise_ =  VectorXd(4,1);
    process_noise_ << 0, 0, 0, 0;

    measurement_noise_ =  VectorXd(4,1);
    measurement_noise_ << 0, 0, 0, 0;

}

KalmanFilter::~KalmanFilter() = default;

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
    // predict the state
    x_ = F_*x_; // + process_noise_;
    P_ = F_*P_*F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
    // update the state by using Kalman Filter equations (LIDAR)
    VectorXd y_;
    y_ = z - H_ * x_;

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * S.inverse();

    //new estimate
    x_ = x_ + (K * y_) + this->process_noise_;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    //update the state by using Extended Kalman Filter equations (RADAR)

    // calculate Jacobian Matrix

    MatrixXd Hj(3,4);

    //recover state parameters
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px * px + py * py;
    double c2 = sqrt(c1);

    //check division by zero
    if(fabs(c1) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    }
    else {
        //compute the Jacobian matrix
        Hj << tools.CalculateJacobian(x_);
    }

    VectorXd y_;
    VectorXd hj(3);
    hj << c2, std::atan2(py, px), (px*vx + py*vy)/c2;
    y_ = z - hj;
    y_(2) = std::fmod(y_(2),M_PI);

    MatrixXd Hjt = Hj.transpose();
    MatrixXd S = Hj * P_ * Hjt + R_;
    MatrixXd PHjt = P_ * Hjt;
    MatrixXd K = PHjt * S.inverse();

    //new estimate
    x_ = x_ + (K * y_);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj) * P_;

}
