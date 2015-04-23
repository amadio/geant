#ifndef SYSTEMOFUNITS_H
#define SYSTEMOFUNITS_H

#include "core/base/Global.h"

// system of units

namespace vecphys {

// length
VECPHYS_GLOBAL double mm = 1.0;
VECPHYS_GLOBAL double cm = 10.*mm;
VECPHYS_GLOBAL double  m = 1.e+3*mm;

// time
VECPHYS_GLOBAL double ns = 1.0;
VECPHYS_GLOBAL double  s = 1.e+9*ns;
VECPHYS_GLOBAL double ms = 1.e-3*s;

// energy
VECPHYS_GLOBAL double MeV = 1.0;
VECPHYS_GLOBAL double TeV = 1.e+6*MeV;
VECPHYS_GLOBAL double GeV = 1.e+3*MeV;
VECPHYS_GLOBAL double KeV = 1.e-3*MeV;
VECPHYS_GLOBAL double  eV = 1.e-6*MeV;

} // end namespace vecphys

#endif
