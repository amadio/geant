#ifndef GEANTV_POSITRONTO2GAMMATESTCOMMON_H
#define GEANTV_POSITRONTO2GAMMATESTCOMMON_H

#include "Common.h"
#include "Geant/PositronTo2GammaModel.h"

using geantphysics::PositronTo2GammaModel;

const double kPos2GammaMinEn = 100 * geant::units::eV;
const double kPos2GammaMaxEn = 100.0 * geant::units::TeV;

#endif
