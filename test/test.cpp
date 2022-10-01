#include <Eigen/Core>
#include <iostream>

int main(int argc, char* argv[]) {
    Eigen::VectorXd rho_e_xn = Eigen::VectorXd::Zero(5);
    Eigen::VectorXd rho_e_xn2 = Eigen::VectorXd::Zero(5);
    for (double i = 5 - 0.75; i >= 0.25; i--) {
        rho_e_xn[5-static_cast<int>(i)-1] = i / 5;
    }

    for (double i = 0.25; i <= 5 - 0.75; i++)
        rho_e_xn2[static_cast<int>(i)] = i / 5;

    std::cout<<rho_e_xn<<std::endl;
    std::cout<<rho_e_xn2<<std::endl;

    return 0;
}