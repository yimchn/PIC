#include <math.h>

#include <iostream>

#include "domain.h"
#include "geometry.h"
#include "output.h"
#include "solver.h"

int main(int argc, char* argv[]) {
    // spdlog::set_level(spdlog::level::debug);
    spdlog::set_level(spdlog::level::info);

    double f = 3e5;  // frenquency, Hz
    double I = 2000 / (Const::PI * pow(0.005, 2));

    // Geometry geo;
    Geometry geo(101, 101, 10, 10, 10, 10);
    geo.SetExtents({-0.5, -0.5}, {0.5, 0.5});

    double dt = geo.dh[0] / (Const::C * sqrt(2));
    int step = static_cast<int>(1 / (f * dt));

    Domain dm(geo, I, f);
    dm.setTime(dt, step);
    Solver solver(dm, 1e5, 1e-6);

    Output<Solver> out(dm, 1e5, 1e-6);

    out.Launch();

    return 0;
}