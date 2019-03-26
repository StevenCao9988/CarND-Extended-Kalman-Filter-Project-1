#include "tools.h"
#include <iostream>

#define ES 0.0001      // 
#define ES2 0.0000001  // 

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.  计算均方根误差
   */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	// check the validity of input  输入数据的合法性检测，增加程序鲁棒性
	// The estimation vector size should not be zero 测量数据的长度不能为0
	if (estimations.size() == 0) {
		cout << "Inout is empty!" << endl;
		return rmse;
	}
	// The estimation vector size should equal ground truth vector size 测量数据长度应该与真实值的长度一致
	if (estimations.size() != ground_truth.size()) {
		cout << "estimations and ground_truth should be the same size!" << endl;
		return rmse;
	}
	// Accumulate squared residuals 累积误差的平方
	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// calculate the mean 计算平均值
	rmse = rmse / estimations.size(); 
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	 * TODO:
	 * Calculate a Jacobian here.  雅克比矩阵
	 */

	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);
	MatrixXd Hj(3, 4);

	// Deal with the special case problems  程序的鲁棒性
	if (fabs(px) < ES and fabs(py) < ES) {
		px = ES;
		py = ES;
	}
	// Pre-compute a set of terms to avoid repeated calculation  提前计算一些复杂的数
	double  QudSum = px * px + py * py;  // QudSum: quadratic sum 平方和
	if (fabs(QudSum) < ES2){     // 程序的鲁棒性
		QudSum = ES2;
	}
	double QudSum_1 = sqrt(QudSum);      // QudSum的1/2 次方
	double QudSum_3 = QudSum * QudSum_1; // QudSum的3/2 次方
	// Compute the Jacobian matrix 正式计算雅可比矩阵
	Hj << px / QudSum_1, py / QudSum_1, 0, 0,
		-(py / QudSum), px / QudSum, 0, 0,
		py * (vx*py - vy * px) / QudSum_3,  px * (vy * px - vx * py) / QudSum_3, px / QudSum_1, py / QudSum_1;
	return Hj; //返回雅克比矩阵
}
