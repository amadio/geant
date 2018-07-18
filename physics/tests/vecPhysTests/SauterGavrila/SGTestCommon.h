#ifndef GEANTV_SAUTERGAVRILATESTCOMMON_H
#define GEANTV_SAUTERGAVRILATESTCOMMON_H

#include "Common.h"
#include "Geant/SauterGavrilaPhotoElectricModel.h"

using geantphysics::SauterGavrilaPhotoElectricModel;

constexpr double kSauGavminEn = 100 * geant::units::eV; // Should be higher then matcut
constexpr double kSauGavmaxEn = 50.0 * geant::units::keV;

#endif
