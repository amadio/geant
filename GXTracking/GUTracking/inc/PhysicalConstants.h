#ifndef PHYSICALCONSTANTS_H
#define PHYSICALCONSTANTS_H

#include "core/base/Global.h"
#include "SystemOfUnits.h"

// HEP Physical Constants

namespace vecphys {

// speed of light
VECPHYS_GLOBAL double c_light = 2.99792458e+8*m/s;
VECPHYS_GLOBAL double c_squared = c_light*c_light;

//atomic equivalent mass unit
VECPHYS_GLOBAL double electron_mass_c2 = 0.510998910 * MeV;
VECPHYS_GLOBAL double inv_electron_mass_c2 = 1.0/electron_mass_c2;
VECPHYS_GLOBAL double amu_c2 = 931.494028*MeV;
VECPHYS_GLOBAL double amu = amu_c2/c_squared;

} // end namespace vecphys

#endif
