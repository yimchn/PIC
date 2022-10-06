#include <math.h>

#include "const.h"
#include "geometry.h"
#include "species.h"

/*constructor*/
Domain::Domain(Geometry geo, double I, double f)
    : geo(geo),
      I(I),
      f(f),
      phi(geo.ni, geo.nj),
      rho(geo.ni, geo.nj),
      E(geo.ni, geo.nj),
      H(geo.ni, geo.nj),
      J(geo.ni, geo.nj) {
    time_start =
        std::chrono::high_resolution_clock::now();  // save starting time point
}

int Domain::getTs() const { return ts; }

double Domain::getTime() const { return time; }

/*returns elapsed wall time in seconds*/
double Domain::getWallTime() {
    auto time_now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time_delta = time_now - time_start;
    return time_delta.count();
}

double Domain::getDt() const { return dt; }

bool Domain::isLastTimeStep() const { return ts == num_ts - 1; }

void Domain::setTime(double dt, int num_ts) {
    this->dt = dt;
    this->num_ts = num_ts;
}

bool Domain::advanceTime() {
    time += dt;
    ts++;
    return ts <= num_ts;
}

Vec2d Domain::XtoL(Vec2d x) {
    Vec2d lc;
    lc[0] = (x[0] - geo.x0[0]) / geo.dh[0];
    lc[1] = (x[1] - geo.x0[1]) / geo.dh[1];
    return lc;
}

/*computes charge density from rho = sum(charge*den)*/
void Domain::computeChargeDensity(std::vector<Species> &species) {
    rho = 0;
    for (Species &sp : species) {
        if (sp.charge == 0) continue;  // don't bother with neutrals
        rho += sp.charge * sp.den;
    }
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double Domain::getPE() {
    double pe = 0;
    for (int i = 0; i < geo.ni; i++)
        for (int j = 0; j < geo.nj; j++) {
            Vec3d efn = E[i][j];  // ef at this node
            double ef2 = efn[0] * efn[0] + efn[1] * efn[1];

            pe += ef2 * geo.node_area[i][j];
        }
    return 0.5 * Const::EPS_0 * pe;
}
