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
   * TODO: Calculate the RMSE here.  ������������
   */
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;
	// check the validity of input  �������ݵĺϷ��Լ�⣬���ӳ���³����
	// The estimation vector size should not be zero �������ݵĳ��Ȳ���Ϊ0
	if (estimations.size() == 0) {
		cout << "Inout is empty!" << endl;
		return rmse;
	}
	// The estimation vector size should equal ground truth vector size �������ݳ���Ӧ������ʵֵ�ĳ���һ��
	if (estimations.size() != ground_truth.size()) {
		cout << "estimations and ground_truth should be the same size!" << endl;
		return rmse;
	}
	// Accumulate squared residuals �ۻ�����ƽ��
	for (unsigned int i = 0; i < estimations.size(); i++) {
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	// calculate the mean ����ƽ��ֵ
	rmse = rmse / estimations.size(); 
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
	 * TODO:
	 * Calculate a Jacobian here.  �ſ˱Ⱦ���
	 */

	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);
	MatrixXd Hj(3, 4);

	// Deal with the special case problems  �����³����
	if (fabs(px) < ES and fabs(py) < ES) {
		px = ES;
		py = ES;
	}
	// Pre-compute a set of terms to avoid repeated calculation  ��ǰ����һЩ���ӵ���
	double  QudSum = px * px + py * py;  // QudSum: quadratic sum ƽ����
	if (fabs(QudSum) < ES2){     // �����³����
		QudSum = ES2;
	}
	double QudSum_1 = sqrt(QudSum);      // QudSum��1/2 �η�
	double QudSum_3 = QudSum * QudSum_1; // QudSum��3/2 �η�
	// Compute the Jacobian matrix ��ʽ�����ſɱȾ���
	Hj << px / QudSum_1, py / QudSum_1, 0, 0,
		-(py / QudSum), px / QudSum, 0, 0,
		py * (vx*py - vy * px) / QudSum_3,  px * (vy * px - vx * py) / QudSum_3, px / QudSum_1, py / QudSum_1;
	return Hj; //�����ſ˱Ⱦ���
}
