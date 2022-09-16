#include <iostream>

#include "const.h"
#include "output.h"
#include "solver.h"
#include "species.h"

using namespace std;
using namespace Const;

int main(int argc, char* argv[]) {
    /*initialize domain*/
    // Domain domain(11, 11);
    // domain.setExtents({-0.1, -0.1}, {0.1, 0.1});
    // domain.setTime(2e-10, 10000);

    // /*set up particle species*/
    // vector<Species> species;
    // species.reserve(2);  // pre-allocate space for two species
    // species.emplace_back(Species("O+", 16 * AMU, QE, domain));
    // species.emplace_back(Species("e-", ME, -1 * QE, domain));

    // /*create particles*/
    // Vec2i np_ions_grid = {21, 21};
    // Vec2i np_eles_grid = {11, 11};
    // species[0].loadParticlesBoxQS(domain.getX0(), domain.getXm(), 1e16,
    //                               np_ions_grid);  // ions
    // species[1].loadParticlesBoxQS(domain.getX0(), domain.getXc(), 1e16,
    //                               np_eles_grid);  // electrons

    // /*initialize potential solver and solve initial potential*/
    // Solver solver(domain, 10000, 1e-4);
    // solver.solve();

    // /*obtain initial electric field*/
    // solver.computeEF();

    // /* main loop*/
    // while (domain.advanceTime()) {
    //     /*move particles*/
    //     for (Species& sp : species) {
    //         sp.advance();
    //         sp.computeNumberDensity();
    //     }

    //     /*compute charge density*/
    //     domain.computeChargeDensity(species);

    //     /*update potential*/
    //     solver.solve();

    //     /*obtain electric field*/
    //     solver.computeEF();

    //     /*screen and file output*/
    //     Output::screenOutput(domain, species);
    //     Output::diagOutput(domain, species);

    //     /*periodically write out results*/
    //     if (domain.getTs() % 100 == 0 || domain.isLastTimeStep())
    //         Output::fields(domain, species);
    // }

    // /* grab starting time*/
    // cout << "Simulation took " << domain.getWallTime() << " seconds\n";
    return 0;  // indicate normal exit    return 0;
}