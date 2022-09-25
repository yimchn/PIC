#pragma once

#include <Eigen/Core>
#include <chrono>
#include <vector>

#include "field.h"
#include "geometry.h"

struct Species;

/*defines the computational dm*/
struct Domain {
    Geometry geo;

    Field<double> phi;  // potential
    Field<double> rho;  // charge density
    Field<Vec3d> E;     // electric field components
    Field<Vec3d> H;     // magnetic field componets
    Field<Vec3d> J;     // current density

    Eigen::ArrayXXd Ex = Eigen::ArrayXXd::Zero(geo.ni, geo.nj + 1);
    Eigen::ArrayXXd Ey = Eigen::ArrayXXd::Zero(geo.ni + 1, geo.nj);
    Eigen::ArrayXXd Ez = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);

    Eigen::ArrayXXd Dx = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Dy = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Dz = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);

    Eigen::ArrayXXd Bx = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd By = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Bz = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);

    Eigen::ArrayXXd Hx = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Hy = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Hz = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);

    Eigen::ArrayXXd Jx = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Jy = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);
    Eigen::ArrayXXd Jz = Eigen::ArrayXXd::Zero(geo.ni, geo.nj);

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
};