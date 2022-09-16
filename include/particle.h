#pragma once

#include <vector>

#include "field.h"

/** Data structures for particle storage **/
struct Particle {
    Vec2d pos;  /*position*/
    Vec2d vel;  /*velocity*/
    double mpw; /*macroparticle weight*/

    Particle(Vec2d x, Vec2d v, double mpw);
};
