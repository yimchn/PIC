#pragma once

#include <Eigen/Core>
#include <chrono>
#include <vector>

#include "field.h"
#include "geometry.h"

struct Species;

using matrix = Eigen::ArrayXXd;

/*defines the computational dm*/
struct Domain {
    Geometry geo;

    double I;
    double f;

    Field<double> phi;  // potential
    Field<double> rho;  // charge density
    Field<Vec3d> E;     // electric field components
    Field<Vec3d> H;     // magnetic field componets
    Field<Vec3d> J;     // current density

    matrix Dx = matrix::Zero(geo.ni, geo.nj);
    matrix Dy = matrix::Zero(geo.ni, geo.nj);
    matrix Dz = matrix::Zero(geo.ni, geo.nj);

    matrix Ex = matrix::Zero(geo.ni, geo.nj);
    matrix Ey = matrix::Zero(geo.ni, geo.nj);
    matrix Ez = matrix::Zero(geo.ni, geo.nj);

    matrix Bx = matrix::Zero(geo.ni, geo.nj);
    matrix By = matrix::Zero(geo.ni, geo.nj);
    matrix Bz = matrix::Zero(geo.ni, geo.nj);

    matrix Hx = matrix::Zero(geo.ni, geo.nj);
    matrix Hy = matrix::Zero(geo.ni, geo.nj);
    matrix Hz = matrix::Zero(geo.ni, geo.nj);

    matrix Jx = matrix::Zero(geo.ni, geo.nj);
    matrix Jy = matrix::Zero(geo.ni, geo.nj);
    matrix Jz = matrix::Zero(geo.ni, geo.nj);

    double dt = 0;    // time step length
    double time = 0;  // physical time
    int ts = -1;      // current time step
    int num_ts = 0;   // number of time steps

    std::chrono::time_point<std::chrono::high_resolution_clock>
        time_start;  // time at simulation start

    /*constructor, allocates memory*/
    Domain(Geometry geo, double I, double f);

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
};