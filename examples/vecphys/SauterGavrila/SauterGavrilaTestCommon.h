#ifndef GEANTV_SAUTERGAVRILATESTCOMMON_H
#define GEANTV_SAUTERGAVRILATESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Geant/SauterGavrilaPhotoElectricModel.h"
#include "Geant/VecSauterGavrilaPhotoElectricModel.h"
#include "Common.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::SauterGavrilaPhotoElectricModel;
using geantphysics::VecSauterGavrilaPhotoElectricModel;

const int kMinBasket = 16;
const int kMaxBasket = 256;
const double minEn   = 100 * geant::units::eV;
const double maxEn   = 100.0 * geant::units::MeV;

SauterGavrilaPhotoElectricModel *PrepareSauterGavrilaModel(bool useAlias)
{
  SauterGavrilaPhotoElectricModel *sauterG = new SauterGavrilaPhotoElectricModel();
  sauterG->SetLowEnergyUsageLimit(minEn);
  sauterG->SetHighEnergyUsageLimit(maxEn);
  sauterG->SetUseSamplingTables(useAlias);
  sauterG->Initialize();
  return sauterG;
}

VecSauterGavrilaPhotoElectricModel *PrepareVecSauterGavrilaModel(bool useAlias)
{
  VecSauterGavrilaPhotoElectricModel *sauterG = new VecSauterGavrilaPhotoElectricModel();
  sauterG->SetLowEnergyUsageLimit(minEn);
  sauterG->SetHighEnergyUsageLimit(maxEn);
  sauterG->SetUseSamplingTables(useAlias);
  sauterG->Initialize();
  return sauterG;
}

void PreparePrimaries(std::vector<LightTrack> &output, int N)
{
  geant::VecRngWrapper rng;
  output.clear();
  for (int i = 0; i < N; ++i) {
    LightTrack gamma;
    double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
    double th  = geant::units::kPi * rng.uniform();
    gamma.SetDirX(sin(phi) * cos(th));
    gamma.SetDirY(cos(phi) * cos(th));
    gamma.SetDirZ(sin(th));
    double eKin = minEn + (maxEn - minEn) * rng.uniform();
    gamma.SetKinE(eKin);
    output.push_back(gamma);
  }
}

void PreparePrimaries(LightTrack_v &output, int N)
{
  geant::VecRngWrapper rng;
  output.SetNtracks(N);
  for (int i = 0; i < N; ++i) {
    double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
    double th  = geant::units::kPi * rng.uniform();
    output.SetDirX(sin(phi) * cos(th), i);
    output.SetDirY(cos(phi) * cos(th), i);
    output.SetDirZ(sin(th), i);
    double eKin = minEn + (maxEn - minEn) * rng.uniform();
    output.SetKinE(eKin, i);
    output.SetTrackIndex(i, i);
  }
}

#endif // GEANTV_SAUTERGAVRILATESTCOMMON_H
