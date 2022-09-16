#pragma once

#include "domain.h"
#include "particle.h"

/*species container*/
struct Species {
    Domain &domain;
    const std::string name; /*species name*/
    const double mass;      /*particle mass in kg*/
    const double charge;    /*particle charge in Coulomb*/

    std::vector<Particle> particles; /*contiguous array for storing particles*/
    Field<double> den;               /*number density*/

    Species(std::string name, double mass, double charge, Domain &domain);

    /*returns the number of simulation particles*/
    size_t getNp();

    /*returns the number of real particles*/
    double getRealCount();

    /*returns the species momentum*/
    Vec2d getMomentum();

    /*returns the species kinetic energy*/
    double getKE();

    /*moves all particles using electric field ef[]*/
    void advance();

    /*compute number density*/
    void computeNumberDensity();

    /*adds a new particle*/
    void addParticle(Vec2d pos, Vec2d vel, double mpwt);

    /*random load of particles in a x1-x2 box representing num_den number
     * density*/
    void loadParticlesSquare(Vec2d x1, Vec2d x2, double num_den, int num_mp);

    /*quiet start load of particles in a x1-x2 box representing num_den number
     * density*/
    void loadParticlesBoxQS(Vec2d x1, Vec2d x2, double num_den, Vec2i num_mp);
};