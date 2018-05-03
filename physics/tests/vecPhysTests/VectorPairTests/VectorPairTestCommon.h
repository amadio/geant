#ifndef GEANTV_VECTORPAIRTESTCOMMON_H
#define GEANTV_VECTORPAIRTESTCOMMON_H

#include "Common.h"
#include "Geant/BetheHeitlerPairModel.h"
#include "Geant/PhysicsParameters.h"
#include "Geant/RelativisticPairModel.h"

using geantphysics::RelativisticPairModel;

using geantphysics::BetheHeitlerPairModel;

constexpr double kBHminEn      = 2.0 * geant::units::kElectronMassC2;
constexpr double kBHmaxEn      = 80.0 * geant::units::GeV;
constexpr double kRelPairMinEn = 80.0 * geant::units::GeV;
constexpr double kRelPairMaxEn = 100.0 * geant::units::TeV;

#endif
