#pragma once

#include <chrono>
#include <vector>

#include "field.h"
#include "geometry.h"

struct Species;

/*defines the computational domain*/
struct Domain {
    Geometry geo;

    Field<double> phi;  // potential
    Field<double> rho;  // charge density
    Field<Vec3d> E;     // electric field components
    Field<Vec3d> H;     // magnetic field componets
    Field<Vec3d> J;     // current density

    double dt = 0;    // time step length
    double time = 0;  // physical time
    int ts = -1;      // current time step
    int num_ts = 0;   // number of time steps

    std::chrono::time_point<std::chrono::high_resolution_clock>
        time_start;  // time at simulation start

    /*constructor, allocates memory*/
    Domain(Geometry geo);

    /*functions for accessing time information*/
    int getTs() const;
    double getTime() const;
    double getWallTime(); /*returns wall time in seconds*/
    double getDt() const;
    bool isLastTimeStep() const;

    /*sets time step and number of time steps*/
    void setTime(double dt, int num_ts);

    /*advances to the next time step, returns true as long as more time steps
     * remain*/
    bool advanceTime();

    /*converts physical position to logical coordinate*/
    Vec2d XtoL(Vec2d x);

    /*computes charge density from rho = sum(charge*den)*/
    void computeChargeDensity(std::vector<Species> &species);

    /*returns the system potential energy*/
    double getPE();

    void UpdateBoundary(double I, double f);
};