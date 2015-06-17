#ifndef SYSTEMOFUNITS_H
#define SYSTEMOFUNITS_H

#include "core/base/Global.h"

// system of units

namespace vecphys {

// length
const double mm = 1.0;
const double cm = 10.*mm;
const double  m = 1.e+3*mm;
const double  meter = 1.e+3*mm;

// time
const double ns = 1.0;
const double  s = 1.e+9*ns;
const double ms = 1.e-3*s;
const double second = s;

// energy
const double MeV = 1.0;
const double TeV = 1.e+6*MeV;
const double GeV = 1.e+3*MeV;
const double keV = 1.e-3*MeV;
const double  eV = 1.e-6*MeV;
const double joule = 6.24150e+12*MeV;

// mass 
const double gram = 1.e-3*joule*s*s/(m*m);

// cross section
const double barn = 1.e-28*m*m;
const double microbarn = 1.e-6*barn;

} // end namespace vecphys

#endif
