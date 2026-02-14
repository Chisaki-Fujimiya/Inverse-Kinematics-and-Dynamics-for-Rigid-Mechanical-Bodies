#include <iostream>
#include <cmath>
#include "Eigen/Dense"
#include <array>

using namespace std;
using namespace Eigen;


MatrixXd jacobian(VectorXd parameter) {
    Matrix<double, 2, 1 > jacobian;
    float theta1 = parameter[3];
    float theta2 = parameter[3] + parameter[4];
    float theta3 = parameter[3] + parameter[4] + parameter[5];
    
    jacobian << 
    parameter[0]*cos(theta1)+parameter[1]*cos(theta2)+parameter[2]*cos(theta3),
    parameter[0]*sin(theta1)+parameter[1]*sin(theta2)+parameter[2]*sin(theta3);

    float err = 1e-3;
    for (int i = 0; i < jacobian.size(); i++) {
        if (jacobian(i) < err) {
            jacobian(i) = 0;
        }
    }
    return jacobian;
}

MatrixXd pseudoinverse(VectorXd parameter) {
    MatrixXd J = jacobian(parameter);
    MatrixXd JT = J.transpose();
    MatrixXd JJT = J * JT;
    MatrixXd inv = JJT.inverse();
    //MatrixXd res = inv * JT;
    /*
    cout << J << endl;
    cout << JT << endl;
    */
    return inv;
}

class Model {
    public:
        VectorXd links;
        VectorXd angle;
        VectorXd target;
};

int main() {
    double PI = 3.141592653589;
    VectorXd parameter(6); parameter << 8.0, 5.0, 3.0, 90.0*PI/180.0, 0.0, 0.0;
    
    cout << pseudoinverse(parameter) << endl;;
    return 1; 
}