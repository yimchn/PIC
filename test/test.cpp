#include <spdlog/spdlog.h>

#include <Eigen/Core>
#include <iostream>
#include <ostream>

// using namespace Eigen;
using vector = Eigen::ArrayXd;
using matrix = Eigen::ArrayXXd;
using Eigen::last;
using Eigen::Map;
using Eigen::seq;

int main(int argc, char* argv[]) {
    vector n = vector::Ones(4);
    matrix m = matrix::Random(4, 4);

    matrix tmp = Map<Eigen::Matrix<double, 4, 4>>(n.data());
    std::cout << tmp << std::endl;

    return 0;
}