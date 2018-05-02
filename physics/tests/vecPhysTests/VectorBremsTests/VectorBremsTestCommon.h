#ifndef GEANTV_VECTORBREMSTESTCOMMON_H
#define GEANTV_VECTORBREMSTESTCOMMON_H

#include "Common.h"
#include "Geant/SeltzerBergerBremsModel.h"

#include "Geant/RelativisticBremsModel.h"

using geantphysics::RelativisticBremsModel;
using geantphysics::SeltzerBergerBremsModel;

const double kSBminEn = 1.0 * geant::units::GeV; // Should be higher then matcut
const double kSBmaxEn = 10.0 * geant::units::GeV;

const double kRelBremsMinEn = 10.0 * geant::units::GeV; // Should be higher then matcut
const double kRelBremsMaxEn = 100.0 * geant::units::TeV;

#endif
