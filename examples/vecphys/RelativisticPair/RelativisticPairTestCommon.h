#ifndef GEANTV_RELATIVISTICPAIRTESTCOMMON_H
#define GEANTV_RELATIVISTICPAIRTESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Common.h"
#include "Geant/RelativisticPairModel.h"
#include "Geant/VecRelativisticPairModel.h"
#include "Geant/MaterialCuts.h"
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::RelativisticPairModel;
using geantphysics::VecRelativisticPairModel;
using geantphysics::Material;
using geantphysics::MaterialCuts;

const int kMinBasket = 16;
const int kMaxBasket = 256;
const double minEn   = 80.0 * geant::units::GeV;
const double maxEn   = 100.0 * geant::units::TeV;

RelativisticPairModel *PrepareRPModel(bool useAlias)
{
  RelativisticPairModel *pairModel = new RelativisticPairModel;
  pairModel->SetLowEnergyUsageLimit(minEn);
  pairModel->SetHighEnergyUsageLimit(maxEn);
  pairModel->SetUseSamplingTables(useAlias);
  pairModel->GetListActiveRegions().resize(1);
  pairModel->GetListActiveRegions()[0]      = true;
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  pairModel->Initialize();
  return pairModel;
}

VecRelativisticPairModel *PrepareVecRPModel(bool useAlias)
{
  VecRelativisticPairModel *pairModel = new VecRelativisticPairModel;
  pairModel->SetLowEnergyUsageLimit(minEn);
  pairModel->SetHighEnergyUsageLimit(maxEn);
  pairModel->SetUseSamplingTables(useAlias);
  pairModel->GetListActiveRegions().resize(1);
  pairModel->GetListActiveRegions()[0] = true;
  pairModel->Initialize();
  return pairModel;
}

void PreparePrimaries(std::vector<LightTrack> &output, int N)
{
  auto vol = vecgeom::GeoManager::Instance().FindLogicalVolume("main_vol");
  assert(vol != nullptr);
  MaterialCuts *cuts = (MaterialCuts *)(vol->GetMaterialCutsPtr());
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
    gamma.SetTrackIndex(i);
    gamma.SetMaterialCutCoupleIndex(cuts->GetIndex());
    output.push_back(gamma);
  }
}

void PreparePrimaries(LightTrack_v &output, int N)
{
  auto vol = vecgeom::GeoManager::Instance().FindLogicalVolume("main_vol");
  assert(vol != nullptr);
  MaterialCuts *cuts = (MaterialCuts *)(vol->GetMaterialCutsPtr());
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
    output.SetMaterialCutCoupleIndex(cuts->GetIndex(), i);
  }
}

void PrepareWorld()
{
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

#endif
