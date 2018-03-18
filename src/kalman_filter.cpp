#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
    x_ << 0.0, 0.0, 0.0, 0.0;

}

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
    // predict the state
    x_ = F_*x_;
    P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    // update the state by using Kalman Filter equations (LIDAR)
    y_ = z - H_ * x_;

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * S.inverse();

    //new estimate
    x_ = x_ + (K * y_);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
    TODO:
      * update the state by using Extended Kalman Filter equations (RADAR)
    */

    // calculate Jacobian Matrix

    MatrixXd Hj(3,4);
    //recover state parameters
    double px = x_(0);
    double py = x_(1);
    double vx = x_(2);
    double vy = x_(3);

    //pre-compute a set of terms to avoid repeated calculation
    double c1 = px*px+py*py;
    double c2 = sqrt(c1);
    double c3 = (c1*c2);

    //check division by zero
    if(fabs(c1) < 0.0001){
        std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    }
    else {
        //compute the Jacobian matrix
        Hj << (px / c2),                     (py / c2),                     0,       0,
              -(py / c1),                    (px / c1),                     0,       0,
              py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;



    }


    y_ = z - H_ * x_;

    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd K = P_ * Ht * S.inverse();

    //new estimate
    x_ = x_ + (K * y_);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;

}
