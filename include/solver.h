#pragma once

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "const.h"
#include "domain.h"

using vector = Eigen::ArrayXd;
using matrix = Eigen::ArrayXXd;
using Eigen::last;
using Eigen::seq;

struct Solver {
    Domain &dm;
    unsigned max_solver_it;  // maximum number of solver iterations
    double tolerance;        // solver tolerance

    // FDTD coefficient
    double dt;
    double delta;
    double eps_r;  // permittivity
    double mu_r;   // permeability
    double sigma;  // conductivity
    double sigma_max;
    double boundary_factor;
    double gradient_k;
    double R0;   // reflection coefficient
    double m;    // grading order
    double eta;  // impedance
    double boundary_width;

    double c1;
    double c2;
    double c3;
    double c4;
    double c5;

    double d1;
    double d2;
    double d3;
    double d4;

    double c1x;
    matrix c2x = matrix::Zero(dm.geo.ni, dm.geo.nj);
    matrix c3x = matrix::Zero(dm.geo.ni, dm.geo.nj);
    vector c4x = vector::Zero(dm.geo.ni);
    vector c5x = vector::Zero(dm.geo.ni);

    double c1y;
    matrix c2y = matrix::Zero(dm.geo.ni, dm.geo.nj);
    matrix c3y = matrix::Zero(dm.geo.ni, dm.geo.nj);
    vector c4y = vector::Zero(dm.geo.nj);
    vector c5y = vector::Zero(dm.geo.nj);

    vector d1z = vector::Zero(dm.geo.ni);
    vector d2z = vector::Zero(dm.geo.ni);
    vector d3z = vector::Zero(dm.geo.nj);
    vector d4z = vector::Zero(dm.geo.nj);

    vector gradientC2 = vector::Zero(dm.geo.npml);
    vector gradientC3 = vector::Zero(dm.geo.npml);
    vector gradientC4 = vector::Zero(dm.geo.npml);
    vector gradientC5 = vector::Zero(dm.geo.npml);

    vector gradientD1 = vector::Zero(dm.geo.npml);
    vector gradientD2 = vector::Zero(dm.geo.npml);
    vector gradientD3 = vector::Zero(dm.geo.npml);
    vector gradientD4 = vector::Zero(dm.geo.npml);

    void Init();

    void UpdatePml();

    void EvaluateFdtd();

    void EvaluateSource(double I, double f, double t);

    /*constructor, sets world*/
    Solver(Domain &domain, int max_it, double tol);

    /*solves potential using Gauss-Seidel*/
    bool solve();

    /*computes electric field = -gradient(phi)*/
    void computeEF();
};
