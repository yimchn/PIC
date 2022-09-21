#include <iostream>

#include "const.h"
#include "output.h"
#include "solver.h"
#include "species.h"

using namespace std;
using namespace Const;

int main(int argc, char* argv[]) {
    /*initialize dm*/
    // Domain dm(11, 11);
    // dm.setExtents({-0.1, -0.1}, {0.1, 0.1});
    // dm.setTime(2e-10, 10000);

    // /*set up particle species*/
    // vector<Species> species;
    // species.reserve(2);  // pre-allocate space for two species
    // species.emplace_back(Species("O+", 16 * AMU, QE, dm));
    // species.emplace_back(Species("e-", ME, -1 * QE, dm));

    // /*create particles*/
    // Vec2i np_ions_grid = {21, 21};
    // Vec2i np_eles_grid = {11, 11};
    // species[0].loadParticlesBoxQS(dm.getX0(), dm.getXm(), 1e16,
    //                               np_ions_grid);  // ions
    // species[1].loadParticlesBoxQS(dm.getX0(), dm.getXc(), 1e16,
    //                               np_eles_grid);  // electrons

    // /*initialize potential solver and solve initial potential*/
    // Solver solver(dm, 10000, 1e-4);
    // solver.solve();

    // /*obtain initial electric field*/
    // solver.computeEF();

    // /* main loop*/
    // while (dm.advanceTime()) {
    //     /*move particles*/
    //     for (Species& sp : species) {
    //         sp.advance();
    //         sp.computeNumberDensity();
    //     }

    //     /*compute charge density*/
    //     dm.computeChargeDensity(species);

    //     /*update potential*/
    //     solver.solve();

    //     /*obtain electric field*/
    //     solver.computeEF();

    //     /*screen and file output*/
    //     Output::screenOutput(dm, species);
    //     Output::diagOutput(dm, species);

    //     /*periodically write out results*/
    //     if (dm.getTs() % 100 == 0 || dm.isLastTimeStep())
    //         Output::fields(dm, species);
    // }

    // /* grab starting time*/
    // cout << "Simulation took " << dm.getWallTime() << " seconds\n";
    return 0;  // indicate normal exit    return 0;
}