#pragma once

#include <random>

/*object for sampling random numbers*/
struct Rnd {
    std::mt19937 mt_gen;                              // random number generator
    std::uniform_real_distribution<double> rnd_dist;  // uniform distribution

    // constructor: set initial random seed and distribution limits
    Rnd();
    double operator()() { return rnd_dist(mt_gen); }
};

extern Rnd rnd;