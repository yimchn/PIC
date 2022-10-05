#include <math.h>

#include <iostream>

#include "domain.h"
#include "geometry.h"
#include "output.h"
#include "solver.h"

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
    // spdlog::set_level(spdlog::level::debug);
    spdlog::set_level(spdlog::level::info);

    double f = 3e5;  // frenquency, Hz
    double I = 1600;

    // Geometry geo;
    Geometry geo(101, 101, 10, 10, 10, 10);
    geo.SetExtents({-0.5, -0.5}, {0.5, 0.5});

    double dt = geo.dh[0] / (Const::C * sqrt(2));
    int step = 10 * static_cast<int>(1 / (f * dt));

    cout << dt << endl;
    cout << step << endl;

    Domain dm(geo);
    dm.setTime(dt, 2 * step);

    Solver solver(dm, 10000, 1e-4);

    std::cout << "Calculating..." << std::endl;
    // solver.UpdateBoundary(dm, I, f);
    // while (dm.ts < 2000) {
    //    solver.UpdateBoundary(dm, I, f);
    //    Output::E(dm);
    //    ++dm.ts;
    //}
    while (dm.advanceTime()) {
        Output::ProgressBar(dm);

        solver.UpdateBoundary(dm, I, f);

        if (dm.ts % 1000 == 0) {
            Output::fields(dm);
        }
    }
    std::cout << "\nCalculation complete" << std::endl;

    return 0;
}