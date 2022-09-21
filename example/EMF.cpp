#include <math.h>

#include <iostream>

#include "domain.h"
#include "geometry.h"
#include "output.h"
#include "solver.h"

using namespace std;

int main(int argc, char *argv[]) {
    spdlog::set_level(spdlog::level::debug);
    double f = 3e5;  // frenquency, Hz
                     //    int step = 150000;
    int step = 400;
    double I = 100;

    Geometry geo;
    geo.SetExtents({-0.1, -0.1}, {0.1, 0.1});

    Domain domain(geo);
    domain.setTime(1 / f / step, step);

    Solver solver(domain, 10000, 1e-4);

    std::cout << "Start calculation" << std::endl;
    while (domain.advanceTime()) {
        // std::cout << "Satrat calculate step: " << dm.ts << "/" << step
        //           << std::endl;
        solver.UpdateBoundary(domain, I, f);
        //        solver.UpdateElectromagnetic();
        //        if (dm.ts % 10000 == 0) Output::fields(dm);
        Output::fields(domain);
    }
    std::cout << "Calculation complete" << std::endl;

    return 0;
}