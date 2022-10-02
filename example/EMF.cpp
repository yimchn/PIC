#include <math.h>

#include <iostream>

#include "domain.h"
#include "geometry.h"
#include "output.h"
#include "solver.h"

using namespace std;

int main(int argc, char *argv[]) {
    // spdlog::set_level(spdlog::level::debug);
    spdlog::set_level(spdlog::level::info);

    double f = 3e5;  // frenquency, Hz
    double I = 1600;

    Geometry geo;
    // Geometry geo(52, 52, 12, 12, 12, 12);
    geo.SetExtents({-0.01, -0.01}, {0.01, 0.01});

    double dt = geo.dh[0] / (Const::C * sqrt(2));
    int step = static_cast<int>(1 / (f * dt));

    Domain dm(geo);
    dm.setTime(dt, step);

    Solver solver(dm, 10000, 1e-4);

    std::cout << "Calculating..." << std::endl;
    while (dm.advanceTime()) {
        Output::ProgressBar(dm.ts, step);
        solver.UpdateBoundary(dm, I, f);
        if (dm.ts % 1000 == 0) {
            Output::fields(dm);
        }
        // Output::fields(dm);
    }
    std::cout << "\nCalculation complete" << std::endl;

    return 0;
}