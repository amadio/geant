#ifndef GEANTV_POSITRONTO2GAMMATESTCOMMON_H
#define GEANTV_POSITRONTO2GAMMATESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Geant/PositronTo2GammaModel.h"
#include "Geant/VecPositronTo2GammaModel.h"
#include "Common.h"
#include <random>

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::PositronTo2GammaModel;
using geantphysics::VecPositronTo2GammaModel;

const int kMinBasket = 16;
const int kMaxBasket = 256;
const double minEn   = 100 * geant::units::eV;
const double maxEn   = 100.0 * geant::units::TeV;

PositronTo2GammaModel *PrepareAnihilModel(bool useAlias)
{
  PositronTo2GammaModel *model = new PositronTo2GammaModel();
  model->SetLowEnergyUsageLimit(minEn);
  model->SetHighEnergyUsageLimit(maxEn);
  model->SetUseSamplingTables(useAlias);
  model->Initialize();
  return model;
}

VecPositronTo2GammaModel *PrepareVecAnihilModel(bool useAlias)
{
  VecPositronTo2GammaModel *model = new VecPositronTo2GammaModel();
  model->SetLowEnergyUsageLimit(minEn);
  model->SetHighEnergyUsageLimit(maxEn);
  model->SetUseSamplingTables(useAlias);
  model->Initialize();
  return model;
}

void PreparePrimaries(std::vector<LightTrack> &output, int N)
{
  std::minstd_rand eng;
  std::uniform_real_distribution<double> rng;
  output.clear();

  for (int i = 0; i < N; ++i) {
    LightTrack positron;
    double phi = geant::units::kTwoPi * rng(eng); // NOT uniform on shpere
    double th  = geant::units::kPi * rng(eng);
    positron.SetDirX(sin(phi) * cos(th));
    positron.SetDirY(cos(phi) * cos(th));
    positron.SetDirZ(sin(th));
    double eKin = minEn + (maxEn - minEn) * rng(eng);
    positron.SetKinE(eKin);
    positron.SetMass(geant::units::kElectronMassC2);
    positron.SetTrackIndex(i);
    output.push_back(positron);
  }
}

void PreparePrimaries(LightTrack_v &output, int N)
{
  std::minstd_rand eng;
  std::uniform_real_distribution<double> rng;
  output.SetNtracks(N);
  for (int i = 0; i < N; ++i) {
    double phi = geant::units::kTwoPi * rng(eng); // NOT uniform on shpere
    double th  = geant::units::kPi * rng(eng);
    output.SetDirX(sin(phi) * cos(th), i);
    output.SetDirY(cos(phi) * cos(th), i);
    output.SetDirZ(sin(th), i);
    double eKin = minEn + (maxEn - minEn) * rng(eng);
    output.SetKinE(eKin, i);
    output.SetMass(geant::units::kElectronMassC2, i);
    output.SetTrackIndex(i, i);
  }
}

#endif
