#ifndef GEANTV_KLEINNISHINATESTCOMMON_H
#define GEANTV_KLEINNISHINATESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"
#include "Geant/TaskData.h"
#include "Geant/KleinNishinaComptonModel.h"
#include "Geant/VecKleinNishinaComptonModel.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::KleinNishinaComptonModel;
using geantphysics::VecKleinNishinaComptonModel;
using geant::TaskData;

const int kMinBasket = 16;
const int kMaxBasket = 256;
const double minEn   = 100 * geant::units::eV;
const double maxEn   = 100.0 * geant::units::TeV;

KleinNishinaComptonModel *PrepareKnishinaModel(bool useAlias)
{
  KleinNishinaComptonModel *kNish = new KleinNishinaComptonModel();
  kNish->SetLowEnergyUsageLimit(minEn);
  kNish->SetHighEnergyUsageLimit(maxEn);
  kNish->SetUseSamplingTables(useAlias);
  kNish->Initialize();
  return kNish;
}

VecKleinNishinaComptonModel *PrepareVecKnishinaModel(bool useAlias)
{
  VecKleinNishinaComptonModel *kNish = new VecKleinNishinaComptonModel();
  kNish->SetLowEnergyUsageLimit(minEn);
  kNish->SetHighEnergyUsageLimit(maxEn);
  kNish->SetUseSamplingTables(useAlias);
  kNish->Initialize();
  return kNish;
}

#endif // GEANTV_KLEINNISHINATESTCOMMON_H
