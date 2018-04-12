#ifndef GEANTV_RELATIVISTICBREMSTESTCOMMON_H
#define GEANTV_RELATIVISTICBREMSTESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Common.h"
#include "Geant/RelativisticBremsModel.h"
#include "Geant/VecRelativisticBremsModel.h"
#include "Geant/MaterialCuts.h"
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::RelativisticBremsModel;
using geantphysics::VecRelativisticBremsModel;
using geantphysics::Material;
using geantphysics::MaterialCuts;

const int kMinBasket = 256;
const int kMaxBasket = 256;
const double minEn   = 10.0 * geant::units::GeV; // Should be higher then matcut
const double maxEn   = 100.0 * geant::units::TeV;

RelativisticBremsModel *PrepareRelBremsModel(bool useAlias)
{
  RelativisticBremsModel *model = new RelativisticBremsModel();
  model->SetLowEnergyUsageLimit(minEn);
  model->SetHighEnergyUsageLimit(maxEn);
  model->SetUseSamplingTables(useAlias);
  model->GetListActiveRegions().resize(1);
  model->GetListActiveRegions()[0]          = true;
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  model->Initialize();
  return model;
}

VecRelativisticBremsModel *PrepareVecRelBremsModel(bool useAlias)
{
  VecRelativisticBremsModel *model = new VecRelativisticBremsModel();
  model->SetLowEnergyUsageLimit(minEn);
  model->SetHighEnergyUsageLimit(maxEn);
  model->SetUseSamplingTables(useAlias);
  model->GetListActiveRegions().resize(1);
  model->GetListActiveRegions()[0] = true;
  model->Initialize();
  return model;
}

void PreparePrimaries(std::vector<LightTrack> &output, int N)
{
  auto vol = vecgeom::GeoManager::Instance().FindLogicalVolume("main_vol");
  assert(vol != nullptr);
  MaterialCuts *cuts = (MaterialCuts *)(vol->GetMaterialCutsPtr());
  geant::VecRngWrapper rng;
  output.clear();
  for (int i = 0; i < N; ++i) {
    LightTrack lepton;
    double phi = geant::units::kTwoPi * rng.uniform(); // NOT uniform on shpere
    double th  = geant::units::kPi * rng.uniform();
    lepton.SetDirX(sin(phi) * cos(th));
    lepton.SetDirY(cos(phi) * cos(th));
    lepton.SetDirZ(sin(th));
    double eKin = minEn + (maxEn - minEn) * rng.uniform();
    lepton.SetKinE(eKin);
    lepton.SetTrackIndex(i);
    lepton.SetMaterialCutCoupleIndex(cuts->GetIndex());
    output.push_back(lepton);
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
