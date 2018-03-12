#ifndef GEANTV_KLEINNISHINATESTCOMMON_H
#define GEANTV_KLEINNISHINATESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"
#include "Geant/TaskData.h"
#include "Geant/KleinNishinaComptonModel.h"

using geantphysics::LightTrack;
using geantphysics::KleinNishinaComptonModel;
using geant::TaskData;

const int kMinBasket = 1 << 4;
const int kMaxBasket = 1 << 8;
const double minEn   = 100 * geant::units::eV;
const double maxEn   = 100.0 * geant::units::TeV;

KleinNishinaComptonModel *PrepareKnishinaModel(bool useAlias)
{
  KleinNishinaComptonModel *kNish = new KleinNishinaComptonModel();
  kNish->SetLowEnergyUsageLimit(minEn);
  kNish->SetHighEnergyUsageLimit(maxEn);
  kNish->SetUseSamplingTables(useAlias);
  kNish->Initialize();
}

#endif // GEANTV_KLEINNISHINATESTCOMMON_H
