#ifndef GEANTV_MOLLERBHABHAESTCOMMON_H
#define GEANTV_MOLLERBHABHAESTCOMMON_H

#include "Common.h"
#include "Geant/MollerBhabhaIonizationModel.h"

using geantphysics::MollerBhabhaIonizationModel;

constexpr double kMollBhminEn = 1.0 * geant::units::GeV; // Should be higher then matcut
constexpr double kMollBhmaxEn = 100.0 * geant::units::TeV;

#endif
