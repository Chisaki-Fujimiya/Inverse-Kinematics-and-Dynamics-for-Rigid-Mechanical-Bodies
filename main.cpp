#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include <array>

using namespace std;
using namespace Eigen;

MatrixXd Jacobian(VectorXd parameter) {
    Matrix<double, 2, 3 > jacobian;
    double theta1 = parameter[3];
    double theta2 = parameter[3] + parameter[4];
    double theta3 = parameter[3] + parameter[4] + parameter[5];
    
    jacobian <<
    //dx/dt
    parameter[0]*sin(theta1)+parameter[1]*sin(theta2)+parameter[2]*sin(theta3), 
    parameter[0]*sin(theta1)+parameter[1]*sin(theta2),
    parameter[2]*sin(theta3),
    //dy/dt
    parameter[0]*cos(theta1)+parameter[1]*cos(theta2)+parameter[2]*cos(theta3),
    parameter[0]*cos(theta1)+parameter[1]*cos(theta2),
    parameter[2]*cos(theta3);

    double err = 1e-3;
    for (int i = 0; i < jacobian.size(); i++) {
        if (jacobian(i) < err) {
            jacobian(i) = 0;
        }
    }
    return jacobian;
}

MatrixXd position(VectorXd parameter) {
    Matrix<double, 2,1> position;
    double theta1 = parameter[3]; double theta2 = parameter[3] + parameter[4];
    double theta3 = parameter[3] + parameter[4] + parameter[5];
    
    position<<
    //x
    parameter[0]*cos(theta1)+parameter[1]*cos(theta2)+parameter[2]*cos(theta3),
    //y
    parameter[0]*sin(theta1)+parameter[1]*sin(theta2)+parameter[2]*sin(theta3);
    return position;
}

MatrixXd pseudoinverse(VectorXd parameter, int iteration, MatrixXd initial_p, MatrixXd error) {
    //initiating iteration values
    MatrixXd Jacobian_iter = Jacobian(parameter); MatrixXd initial_iter = initial_p;
    MatrixXd JacobianT = Jacobian_iter.transpose();
    Matrix<double, 3, 1> angle_iter;
    angle_iter << parameter[3], parameter[4], parameter[5];

    //Damping values to avoid singularity
    double alpha = 0.01; double lambda = 1e-5;
    MatrixXd Ide = MatrixXd::Identity(Jacobian_iter.cols(), Jacobian_iter.cols());
    MatrixXd dls = lambda * lambda * Ide;

    //Looping angle to approximate target location with the angle
    MatrixXd res;
    for (int i = 0; i < iteration+1; i ++) {
        VectorXd temp(6); temp << 8.0e+0, 5.0e+0, 3.0e+0, angle_iter[0], angle_iter[1], angle_iter[2];

        res = alpha*((((JacobianT * Jacobian_iter) + dls).inverse() * JacobianT)*error);
        angle_iter = angle_iter + rem;

        Jacobian_iter = Jacobian(temp); 
        JacobianT = Jacobian_iter.transpose();
        initial_iter = position(temp);

        Vector3d angle_deg = angle_iter * 180.0 / M_PI;
        Vector2d current_position = position(temp);
        cout << "Iteration: " << i 
            << " , angles: " << angle_deg.transpose() 
            << " , position: " << current_position.transpose() << endl;
    }
    return res;
}

int main() {
    double PI = 3.141592653589;
    VectorXd parameter(6); parameter << 8.0, 5.0, 3.0, 45.0*PI/180.0, 45.0*PI/180.0, 45.0*PI/180.0;
    MatrixXd initial_p = position(parameter);
    Matrix<double, 2, 1> Target; Target << 10.0, 12.0;
    pseudoinverse(parameter, 100, initial_p, Target - initial_p);
    return 1; 
}

