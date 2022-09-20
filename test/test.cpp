//#include <math.h>

#include <Eigen/Core>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <unsupported/Eigen/CXX11/Tensor>

#include "const.h"
#include "field.h"
#include "utlis.h"

//using namespace Eigen;
using namespace std;

int main(int argc, char* argv[]) {
    // Eigen::ArrayXXd m;
    // Eigen::ArrayXXd n;
    // m.resize(2, 2);
    // n.resize(2, 2);

    // m << 1, 2, 3, 4;
    // n << 1, 2, 3, 4;

    // MatrixX<Vector3d> mat;
    // mat.resize(2, 2);
    // mat(0, 0).setConstant(1);
    // mat(0, 0)(0) = 1;
    // mat(0, 0)(1) = 2;
    // mat(0, 0)(2) = 3;
    // mat(0, 1).setConstant(2);
    // mat(1, 0).setConstant(3);
    // mat(1, 1).setConstant(4);
    // cout << mat(0, 0)(2) << endl;

    //    Tensor<Vector3d, 3> t;
    //    t.resize(2, 2, 2);
    //    t(0, 0, 0)(0) = 1;
    //    t(0, 0, 0)(1) = 2;
    //    t(0, 0, 0)(2) = 3;
    //    std::cout << t(0, 0, 0)(0) << std::endl;
    //    std::cout << t(0, 0, 0)(1) << std::endl;
    //    std::cout << t(0, 0, 0)(2) << std::endl;

    Eigen::ArrayXd m;
    m.resize(4);
    m << 1, 2, 3, 4;
    std::cout << 2*Eigen::pow(m, 2) << std::endl;

    return 0;
}