#ifndef GEANTV_VECTORPAIRTESTCOMMON_H
#define GEANTV_VECTORPAIRTESTCOMMON_H

#include "Common.h"
#include "Geant/BetheHeitlerPairModel.h"
#include "Geant/PhysicsParameters.h"
#include "Geant/RelativisticPairModel.h"

using geantphysics::RelativisticPairModel;

using geantphysics::BetheHeitlerPairModel;

const double kBHminEn      = 2.0 * geant::units::kElectronMassC2;
const double kBHmaxEn      = 80.0 * geant::units::GeV;
const double kRelPairMinEn = 80.0 * geant::units::GeV;
const double kRelPairMaxEn = 100.0 * geant::units::TeV;

#endif
