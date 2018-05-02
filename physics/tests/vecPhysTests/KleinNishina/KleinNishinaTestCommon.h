#ifndef GEANTV_KLEINNISHINATESTCOMMON_H
#define GEANTV_KLEINNISHINATESTCOMMON_H

#include "Common.h"
#include "Geant/KleinNishinaComptonModel.h"

using geantphysics::KleinNishinaComptonModel;

const double kKNminEn = 100 * geant::units::eV;
const double kKNmaxEn = 100.0 * geant::units::TeV;

#endif // GEANTV_KLEINNISHINATESTCOMMON_H
