#ifndef GEANTV_MOLLERBHABHAESTCOMMON_H
#define GEANTV_MOLLERBHABHAESTCOMMON_H

#include "Common.h"
#include "Geant/MollerBhabhaIonizationModel.h"

using geantphysics::MollerBhabhaIonizationModel;

const double kMollBhminEn = 1.0 * geant::units::GeV; // Should be higher then matcut
const double kMollBhmaxEn = 100.0 * geant::units::TeV;

#endif
