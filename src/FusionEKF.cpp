#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices

    //measurement covariance matrix - laser
    R_laser_ = MatrixXd(2, 2);
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    H_laser_ = MatrixXd(2, 4);
    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;


    /**
    TODO:
      * Finish initializing the FusionEKF.
      * Set the process and measurement noises
    */
    ax = ay = 3.0; // maximum acceleration. User-defined. [m/s^2]
    ekf_.process_noise_ << 0, 0, 0, 0;

    ekf_.measurement_noise_ << 0, 0, 0, 0;




}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << std::endl;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */

            double px = measurement_pack.raw_measurements_(0) * std::cos(measurement_pack.raw_measurements_(2));
            double py = measurement_pack.raw_measurements_(0) * std::sin(measurement_pack.raw_measurements_(2));
            ekf_.x_ << px, py, 0, 0;
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */

            double px = measurement_pack.raw_measurements_(0);
            double py = measurement_pack.raw_measurements_(1);
            ekf_.x_ << px, py, 0, 0;
        }

        std::cout << "Initialized!!!"<< std::endl;
        previous_timestamp_ = measurement_pack.timestamp_;
        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix Q.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */
    double delta_t = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
    ekf_.F_(0,2) = delta_t;
    ekf_.F_(1,3) = delta_t;
    //ekf_.F_ << 1, 0, delta_t, 0,
    //           0, 1, 0,       delta_t,
    //           0, 0, 1,       0,
    //           0, 0, 0,       1;

    double t4_ = std::pow(delta_t, 4) / 4.0;
    double t3_ = std::pow(delta_t, 3) / 2.0;
    double t2_ = std::pow(delta_t, 2);
    double noise_ax2 = 9.0*9.0;
    double noise_ay2 = 9.0*9.0;
    ekf_.Q_ << t4_* noise_ax2, 0, t3_ * noise_ax2, 0,
            0, t4_*noise_ay2, 0, t3_*noise_ay2,
            t3_*noise_ax2, 0, t2_*noise_ax2, 0,
            0, t3_*noise_ay2, 0, t2_*noise_ay2;

    ekf_.process_noise_ << (ax*t2_)/2, (ay*t2_)/2, ax*delta_t, ay*delta_t;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
