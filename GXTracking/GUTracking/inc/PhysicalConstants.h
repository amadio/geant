#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

#include "core/base/Global.h"
#include "SystemOfUnits.h"

// HEP Physical Constants

namespace vecphys {

// speed of light
const double c_light = 2.99792458e+8*m/s;
const double c_squared = c_light*c_light;

//atomic equivalent mass unit
const double electron_mass_c2 = 0.510998910 * MeV;
const double inv_electron_mass_c2 = 1.0/electron_mass_c2;
const double amu_c2 = 931.494028*MeV;
const double amu = amu_c2/c_squared;

} // end namespace vecphys

#endif
