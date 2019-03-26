#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"
#define EPS 0.0001 // A very small number

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

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
		double rho = measurement_pack.raw_measurements_[0];         // range 轴半径
		double phi = measurement_pack.raw_measurements_[1];         // bearing 偏转角 
		double rho_dot = measurement_pack.raw_measurements_[2];     // velocity of rho 轴向速度 
      //Convert radar from polar to cartesian coordinates  极坐标转换成直角坐标
		double x = rho * cos(phi);
		double y = rho * sin(phi);
		double vx = rho_dot * cos(phi);
		double vy = rho_dot * sin(phi);
		ekf_.x_ << x, y, vx, vy;
    }else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
	  // 激光雷达只能测量位置，并不直接测量速度,所以给0
		ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
	// Deal with the special case initialisation problems 鲁棒性处理
	if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS) {
		ekf_.x_(0) = EPS;
		ekf_.x_(1) = EPS;
	}
	// Initial covariance matrix 初始化协方差矩阵
	ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ << 1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1000, 0,
		0, 0, 0, 1000;
	// Print the initialization results 打印卡尔曼滤波初始数据
	cout << "EKF init: " << ekf_.x_ << endl;
	// Save the initiall timestamp for dt calculation 记录每一次的时间 用于预测和迭代
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
  // Calculate the timestep between measurements in seconds 计算两次观测之间的时间差
  float dt = (measurement_pack.timestamp_ - previous_timestamp_);
  if (dt > 0){
	  dt /= 1000000.0;  // 微妙转化成秒
  }else{
	  cout << "error time data!" << endl;
  }

  previous_timestamp_ = measurement_pack.timestamp_;
  // State transition matrix update 状态转换矩阵更新
  ekf_.F_ = MatrixXd(4, 4);
      // v * dt = dl 速度乘以时间 等于距离
  ekf_.F_ << 1, 0, dt, 0,
	  0, 1, 0, dt,
	  0, 0, 1, 0,
	  0, 0, 0, 1;

  // Noise covariance matrix computation  噪音协方差矩阵计算
  // Noise values from the task 根据任务确定噪音值
  float noise_ax = 9.0;
  float noise_ay = 9.0;
  // Precompute some usefull values to speed up calculations of Q 提前计算一些对矩阵Q有用的量
  double dt_2 = dt * dt;   // dt^2
  double dt_3 = dt_2 * dt; // dt^3
  double dt_4 = dt_3 * dt; // dt^4
  double dt_4_4 = dt_4 / 4; // dt^4 / 4
  double dt_3_2 = dt_3 / 2; // dt^3 / 2
  ekf_.Q_ = MatrixXd(4, 4); 
  // 以加速度为基准的动力学方程 计算距离与速度之间的协方差矩阵
  ekf_.Q_ << dt_4_4 * noise_ax, 0, dt_3_2 * noise_ax, 0,
	  0, dt_4_4 * noise_ay, 0, dt_3_2 * noise_ay,
	  dt_3_2 * noise_ax, 0, dt_2 * noise_ax, 0,
	  0, dt_3_2 * noise_ay, 0, dt_2 * noise_ay;

  // State transition 
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
	  // 雅克比矩阵	  
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // TODO: Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
	  //cout << "x_ = " << ekf_.x_ << endl;
  } else {
	  cout << "error data" << endl;
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
