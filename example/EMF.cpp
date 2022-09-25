#include <math.h>

#include <iostream>

#include "domain.h"
#include "geometry.h"
#include "output.h"
#include "solver.h"

using namespace std;

int main(int argc, char *argv[]) {
    spdlog::set_level(spdlog::level::info);
    //    spdlog::set_level(spdlog::level::debug);
    double f = 3e5;  // frenquency, Hz
    int step = 150000;
    // int step = 400;
    double I = 100;

    Geometry geo;
    geo.SetExtents({-0.1, -0.1}, {0.1, 0.1});

    Domain dm(geo);
    dm.setTime(1 / f / step, step);

    Solver solver(dm, 10000, 1e-4);

    std::cout << "Start calculation" << std::endl;
    while (dm.advanceTime()) {
        spdlog::info("Calculating: {} / {}", dm.ts, dm.num_ts);
        solver.EvaluateSource(I, f, dm.getTime());
        solver.EvaluateFdtd();
        if (dm.ts % 1000 == 0) Output::fields(dm);
    }
    std::cout << "Calculation complete" << std::endl;

    return 0;
}