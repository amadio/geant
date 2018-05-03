#ifndef GEANTV_VECTORBREMSTESTCOMMON_H
#define GEANTV_VECTORBREMSTESTCOMMON_H

#include "Common.h"
#include "Geant/SeltzerBergerBremsModel.h"

#include "Geant/RelativisticBremsModel.h"

using geantphysics::RelativisticBremsModel;
using geantphysics::SeltzerBergerBremsModel;

constexpr double kSBminEn = 1.0 * geant::units::GeV; // Should be higher then matcut
constexpr double kSBmaxEn = 10.0 * geant::units::GeV;

constexpr double kRelBremsMinEn = 10.0 * geant::units::GeV; // Should be higher then matcut
constexpr double kRelBremsMaxEn = 100.0 * geant::units::TeV;

#endif
