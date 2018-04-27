#ifndef GEANTV_COMMON_H
#define GEANTV_COMMON_H

#include <Geant/Gamma.h>
#include <Geant/Electron.h>
#include <Geant/Positron.h>
#include "Geant/PhysicsData.h"
#include "Geant/LightTrack.h"

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Common.h"
#include "Geant/MaterialCuts.h"
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"

#include <vector>
#include <array>
#include <Geant/EMModel.h>

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::EMModel;
using geantphysics::Material;
using geantphysics::MaterialCuts;
using geant::TaskData;

const int kBasketSize    = 256;
const int kBasketSamples = 10000;

std::unique_ptr<TaskData> PrepareTaskData()
{
  auto PhysData    = new geantphysics::PhysicsData();
  auto Td          = std::unique_ptr<TaskData>(new TaskData(1, 2048));
  Td->fPhysicsData = PhysData;
  return Td;
}

void CleanTaskData(TaskData *td)
{
  if (td->fPhysicsData) delete td->fPhysicsData;
}

const double kWrongVal = 1000.0;

double XYDirToAzumuth(double Xdir, double Ydir)
{
  double az = kWrongVal;
  if (Xdir != 0.0 || Ydir != 0.0) {
    az = std::atan2(Xdir, Ydir) / 3.1416;
  }
  return az;
}

double ZDirLog(double ZDir)
{
  double logElz = 0.5 * (1.0 - ZDir);
  if (logElz > 0 && log10(logElz) > -12.0) {
    logElz = log10(logElz);
  } else {
    logElz = kWrongVal;
  }
  return logElz;
}

void PrepareWorld()
{
  // TODO

  Printf("Preparing world");
  vecgeom::GeoManager::Instance().Clear();
  std::string materialName("NIST_MAT_Pb"); // material is lead
  Material *matDetector = Material::NISTMaterial(materialName);

  // Set production cuts if needed: not needed in case of pair-production
  bool iscutinlength  = true;
  double prodCutValue = 1. * geant::units::mm;
  double gcut         = prodCutValue;
  double emcut        = prodCutValue;
  double epcut        = prodCutValue;

  vecgeom::UnplacedBox worldParams = vecgeom::UnplacedBox(1., 1., 1.);
  vecgeom::LogicalVolume *worldl   = new vecgeom::LogicalVolume("main_vol", &worldParams);
  // create one region and assigne to the logical volume
  vecgeom::Region *aRegion = new vecgeom::Region("ARegion", iscutinlength, gcut, emcut, epcut);
  worldl->SetRegion(aRegion);
  // set the material pointer in the world logical volume
  worldl->SetMaterialPtr((void *)matDetector);
  vecgeom::GeoManager::Instance().SetWorld(worldl->Place());
  vecgeom::GeoManager::Instance().CloseGeometry();
  // Create all(we have only one) MaterialCuts
  geantphysics::MaterialCuts::CreateAll();
}

enum class TestParticleType {
  Gamma,
  Em,
  Ep,
};

void CreateParticles(double minEn, double maxEn, bool randomDir, TestParticleType type, LightTrack_v &out, int N)
{
  geant::RngWrapper rng;
  out.SetNtracks(N);
  for (int i = 0; i < N; ++i) {
    double E = 0.0;
    if (minEn == maxEn) {
      E = minEn;
    } else {
      E = rng.uniform(minEn, maxEn);
    }
    out.SetKinE(E, i);

    if (randomDir) {
      double dX = rng.Gauss(0.0, 1.0);
      double dY = rng.Gauss(0.0, 1.0);
      double dZ = rng.Gauss(0.0, 1.0);
      double dL = std::sqrt(dX * dX + dY * dY + dZ * dZ);
      out.SetDirX(dX / dL, i);
      out.SetDirY(dY / dL, i);
      out.SetDirZ(dZ / dL, i);
    } else {
      out.SetDirX(0.0, i);
      out.SetDirY(0.0, i);
      out.SetDirZ(1.0, i);
    }

    geantphysics::Particle *particle;
    switch (type) {
    case TestParticleType::Gamma:
      particle = geantphysics::Gamma::Definition();
      break;
    case TestParticleType::Em:
      particle = geantphysics::Electron::Definition();
      break;
    case TestParticleType::Ep:
      particle = geantphysics::Positron::Definition();
      break;
    }
    out.SetMass(particle->GetPDGMass(), i);
    out.SetGVcode(particle->GetInternalCode(), i);

    out.SetEnergyDeposit(0.0, i);
    out.SetTrackIndex(i, i);

    std::vector<vecgeom::LogicalVolume *> logicalVols;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(logicalVols);
    assert(logicalVols.size() != 0);
    int volIdx         = (int)rng.uniform() * logicalVols.size();
    auto vol           = logicalVols[volIdx];
    MaterialCuts *cuts = (MaterialCuts *)(vol->GetMaterialCutsPtr());
    out.SetMaterialCutCoupleIndex(cuts->GetIndex(), i);
  }
}

std::vector<double> GetE(const LightTrack_v &tracks, bool total = false)
{
  std::vector<double> energy;
  energy.resize(tracks.GetNtracks());
  for (int i = 0; i < tracks.GetNtracks(); ++i) {
    energy[i] = tracks.GetKinE(i);
    if (total) {
      energy[i] += tracks.GetEnergyDeposit(i) + tracks.GetMass(i);
    }
  }
  return energy;
}

std::vector<double> GetE(const std::vector<LightTrack> tracks, bool total = false)
{
  std::vector<double> energy;
  energy.resize(tracks.size());
  for (size_t i = 0; i < tracks.size(); ++i) {
    energy[i] = tracks[i].GetKinE();
    if (total) {
      energy[i] += tracks[i].GetEnergyDeposit() + tracks[i].GetMass();
    }
  }
  return energy;
}

double GetTotalE(const LightTrack_v &tracks)
{
  double E = 0.0;
  for (int i = 0; i < tracks.GetNtracks(); ++i) {
    E += tracks.GetKinE(i) + tracks.GetMass(i) + tracks.GetEnergyDeposit(i);
  }
  return E;
}

double GetTotalE(const LightTrack *tracks, int N)
{
  double E = 0.0;
  for (int i = 0; i < N; ++i) {
    E += tracks[i].GetKinE() + tracks[i].GetMass() + tracks[i].GetEnergyDeposit();
  }
  return E;
}

std::array<double, 3> GetTotalP(const LightTrack_v &tracks)
{
  std::array<double, 3> Pvec = {0.0, 0.0, 0.0};
  for (int i = 0; i < tracks.GetNtracks(); ++i) {
    double T = tracks.GetKinE(i);
    double m = tracks.GetMass(i);
    double P = std::sqrt(T * (T + 2.0 * m));
    Pvec[0] += P * tracks.GetDirX(i);
    Pvec[1] += P * tracks.GetDirY(i);
    Pvec[2] += P * tracks.GetDirZ(i);
  }
  return Pvec;
}

std::array<double, 3> GetTotalP(const LightTrack *tracks, int N)
{
  std::array<double, 3> Pvec = {0.0, 0.0, 0.0};
  for (int i = 0; i < N; ++i) {
    double T = tracks[i].GetKinE();
    double m = tracks[i].GetMass();
    double P = std::sqrt(T * (T + 2.0 * m));
    Pvec[0] += P * tracks[i].GetDirX();
    Pvec[1] += P * tracks[i].GetDirY();
    Pvec[2] += P * tracks[i].GetDirZ();
  }
  return Pvec;
}

void CheckEnergyMomentumConservation(double E0, std::array<double, 3> P0, const LightTrack_v &scatteredPrims,
                                     const LightTrack_v &secondaries)
{
  double Eprim1 = GetTotalE(scatteredPrims);
  auto Pprim1   = GetTotalP(scatteredPrims);
  double Esec1  = GetTotalE(secondaries);
  auto Psec1    = GetTotalP(secondaries);
  if (std::abs(E0 - (Eprim1 + Esec1)) / E0 > 1.e-5) {
    Printf("Energy is not conserved: %f %f", E0, Eprim1 + Esec1);
  }
  double P0Abs = P0[0] * P0[0] + P0[1] * P0[1] + P0[2] * P0[2];
  for (int i = 0; i < 3; ++i) {
    if (std::abs(P0[i] - (Pprim1[i] + Psec1[i])) / P0Abs > 1.e-5) {
      Printf("Momentum is not conserved: %f %f", P0[i], Pprim1[i] + Psec1[i]);
    }
  }
}

void CheckEnergyMomentumConservation(double E0, std::array<double, 3> P0, const LightTrack *scatteredPrims, int primN,
                                     const LightTrack *secondaries, int secN)
{
  double Eprim1 = GetTotalE(scatteredPrims, primN);
  auto Pprim1   = GetTotalP(scatteredPrims, primN);
  double Esec1  = GetTotalE(secondaries, secN);
  auto Psec1    = GetTotalP(secondaries, secN);
  if (std::abs(E0 - (Eprim1 + Esec1)) / E0 > 1.e-5) {
    Printf("Energy is not conserved: %f %f", E0, Eprim1 + Esec1);
  }
  double P0Abs = P0[0] * P0[0] + P0[1] * P0[1] + P0[2] * P0[2];
  for (int i = 0; i < 3; ++i) {
    if (std::abs(P0[i] - (Pprim1[i] + Psec1[i])) / P0Abs > 1.e-5) {
      Printf("Momentum is not conserved: %f %f", P0[i], Pprim1[i] + Psec1[i]);
    }
  }
}

std::vector<int> GetGVCodeHist(const LightTrack_v &tracks)
{
  std::vector<int> ptHist;
  auto &ptTable = geantphysics::Particle::GetInternalParticleTable();
  ptHist.insert(ptHist.begin(), ptTable.size(), 0);
  for (int i = 0; i < tracks.GetNtracks(); ++i) {
    assert(tracks.GetGVcode(i) >= 0 && tracks.GetGVcode(i) < ptTable.size());
    ptHist[tracks.GetGVcode(i)]++;
  }
  return ptHist;
}

std::vector<int> GetGVCodeHist(const LightTrack *tracks, int N)
{
  std::vector<int> ptHist;
  auto &ptTable = geantphysics::Particle::GetInternalParticleTable();
  ptHist.insert(ptHist.begin(), ptTable.size(), 0);
  for (int i = 0; i < N; ++i) {
    assert(tracks[i].GetGVcode() >= 0 && tracks[i].GetGVcode() < ptTable.size());
    ptHist[tracks[i].GetGVcode()]++;
  }
  return ptHist;
}

void CreateParticles(double minEn, double maxEn, bool randomDir, TestParticleType type, std::vector<LightTrack> &out,
                     int N)
{
  out.resize(N);
  geant::RngWrapper rng;
  for (int i = 0; i < N; ++i) {
    double E = 0.0;
    if (minEn == maxEn) {
      E = minEn;
    } else {
      E = rng.uniform(minEn, maxEn);
    }
    out[i].SetKinE(E);

    if (randomDir) {
      double dX = rng.Gauss(0.0, 1.0);
      double dY = rng.Gauss(0.0, 1.0);
      double dZ = rng.Gauss(0.0, 1.0);
      double dL = std::sqrt(dX * dX + dY * dY + dZ * dZ);
      out[i].SetDirX(dX / dL);
      out[i].SetDirY(dY / dL);
      out[i].SetDirZ(dZ / dL);
    } else {
      out[i].SetDirX(0.0);
      out[i].SetDirY(0.0);
      out[i].SetDirZ(1.0);
    }

    geantphysics::Particle *particle;
    switch (type) {
    case TestParticleType::Gamma:
      particle = geantphysics::Gamma::Definition();
      break;
    case TestParticleType::Em:
      particle = geantphysics::Electron::Definition();
      break;
    case TestParticleType::Ep:
      particle = geantphysics::Positron::Definition();
      break;
    }
    out[i].SetMass(particle->GetPDGMass());
    out[i].SetGVcode(particle->GetInternalCode());

    out[i].SetEnergyDeposit(0.0);
    out[i].SetTrackIndex(i);

    std::vector<vecgeom::LogicalVolume *> logicalVols;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(logicalVols);
    assert(logicalVols.size() != 0);
    int volIdx         = (int)rng.uniform() * logicalVols.size();
    auto vol           = logicalVols[volIdx];
    MaterialCuts *cuts = (MaterialCuts *)(vol->GetMaterialCutsPtr());
    out[i].SetMaterialCutCoupleIndex(cuts->GetIndex());
  }
}

void SampleSecondariesScalar(geantphysics::EMModel *model, std::vector<LightTrack> &primaries, geant::TaskData *td)
{
  td->fPhysicsData->ClearSecondaries();
  for (size_t i = 0; i < primaries.size(); ++i) {
    model->SampleSecondaries(primaries[i], td);
  }
}

void SampleSecondariesVector(geantphysics::EMModel *model, LightTrack_v &primaries, geant::TaskData *td)
{
  td->fPhysicsData->GetSecondarySOA().ClearTracks();
  model->SampleSecondariesVector(primaries, td);
}

std::unique_ptr<EMModel> InitEMModel(geantphysics::EMModel *model, double minEn, double maxEn, double alias)
{
  model->SetLowEnergyUsageLimit(minEn);
  model->SetHighEnergyUsageLimit(maxEn);
  model->SetUseSamplingTables(alias);
  model->GetListActiveRegions().resize(1);
  model->GetListActiveRegions()[0]          = true;
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  model->Initialize();
  return std::unique_ptr<EMModel>(model);
}

#endif // GEANTV_COMMON_H
