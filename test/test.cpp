#include <math.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>

#include "const.h"
#include "field.h"
#include "utlis.h"

void change(int (&a)[5]);

int main(int argc, char *argv[]) {
    // output file
    // std::ofstream out("../results/hello.txt");
    // if (out.is_open()) {
    //     out << "This is a line.\n";
    //     out << "This is another line.\n";
    //     out.close();
    //     std::cout << "ok" << std::endl;
    // } else {
    //     std::cout << "failed" << std::endl;
    // }

    // field
    // double f = 3e5;
    // double step = 400;
    // double I = 1600;
    // double dt = 1 / f / step;
    // int s = 4;
    // Rnd rnd;
    // Field<Vec3d> ef(s, s);

    // for (int i = 0; i < s; i++) {
    //     for (int j = 0; j < s; j++) {
    //         // ef[i][j] = rnd();
    //         ef[i][j] = 0;
    //     }
    // }

    // for (int ns = 0; ns < step; ns++) {
    //     dt += dt;
    //     for (int i = 0; i < s; i++) {
    //         for (int j = 0; j < s; j++) {
    //             ef[1][j] = I * sin(2 * Const::PI * f * dt);
    //             ef[s - 2][j] = I * sin(2 * Const::PI * f * dt);
    //             ef[i][1] = I * cos(2 * Const::PI * f * dt);
    //             ef[i][s - 2] = I * cos(2 * Const::PI * f * dt);

    //             if (j == s - 1)
    //                 std::cout << std::setprecision(2) << ef[i][j] << "\n";
    //             else
    //                 std::cout << std::setprecision(2) << ef[i][j] << "\t";
    //         }
    //     }
    //     std::cout << "------------------------------------------------------"
    //               << std::endl;
    // }

    // std::cout << 1 / 3e5 / 400 / (Const::EPS_0 * 0.02) << std::endl;
    // std::cout << 1 / 3e5 / 400 << std::endl;

    double n = 20;
    n -= 0.75;
    std::cout << static_cast<int>(n) << std::endl;

    return 0;
}

void change(int (&a)[5]) {
    for (int i = 0; i < 5; i++) {
        a[i] = i + 1;
    }
}