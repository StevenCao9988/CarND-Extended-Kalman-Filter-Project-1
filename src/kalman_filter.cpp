#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

#define EPS 0.001;  // A very small number
#define EPS2 0.000001; // A very small number;
#define YMIN 0.1;    // y ����Сֵ

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
  // ��ͨ���ԵĿ������˲� ���϶����ֳɵ�
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;  // ������ϵ��

	if (x_(0) < 0)
	{
		;
	}
	x_ = x_ + (K * y);      
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
	// Recalculate x object state to rho, phi, rho_dot coordinates ֱ�����껹ԭ�ɼ�����
	double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));           // range ��뾶
	double phi = atan2(x_(1), x_(0));                       // bearing ƫת��  ������atan2�������޷�����-�� ~ +��
	double rho_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;     // velocity of rho �����ٶ� 
	VectorXd h = VectorXd(3);
	h << rho, phi, rho_dot;

	VectorXd y = z - h;
	
    if (y(1) > 0.1)
	{
		y(1) = 0.1;
	}
	else if (y(1) < (-0.1))
	{
	    y(1) = (-0.1);
	}
	cout << "z_ = " << y << endl;
	cout << "h_ = " << h << endl;
	cout << "y_ = " << y << endl;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;  // ������ϵ��

	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
