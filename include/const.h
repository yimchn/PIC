#pragma once

#include <math.h>

namespace Const {
const double EPS_0 = 8.85418782e-12;    // C/(V*m), vacuum permittivity
const double QE = 1.602176565e-19;      // C, electron charge
const double AMU = 1.660538921e-27;     // kg, atomic mass unit
const double ME = 9.10938215e-31;       // kg, electron mass
const double K = 1.380648e-23;          // J/K, Boltzmann constant
const double PI = acos(-1);             // pi
const double EvToK = QE / K;            // 1eV in K ~ 11604
const double C = 2.99792458e8;          // speed of light
const double MU_0 = 4.0 * PI * 1.0e-7;  // H/m, Vacuum permeability
}  // namespace Const
