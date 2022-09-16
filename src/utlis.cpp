#include "utlis.h"

 Rnd::Rnd() : mt_gen{std::random_device()()}, rnd_dist{0, 1.0} {}
