#ifndef GEANTV_BETHEHEITLERTESTCOMMON_H
#define GEANTV_BETHEHEITLERTESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Common.h"
#include "Geant/BetheHeitlerPairModel.h"
#include "Geant/VecBetheHeitlerPairModel.h"
#include "Geant/MaterialCuts.h"
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::BetheHeitlerPairModel;
using geantphysics::VecBetheHeitlerPairModel;
using geantphysics::Material;
using geantphysics::MaterialCuts;

const int kMinBasket      = 16;
const int kMaxBasket      = 256;
const double minEn        = 100.0 * geant::units::eV;
const double minPrimaryEn = 2.0 * geant::units::kElectronMassC2;
const double maxEn        = 80.0 * geant::units::GeV;

BetheHeitlerPairModel *PrepareBHModel(bool useAlias)
{
  BetheHeitlerPairModel *BHmodel = new BetheHeitlerPairModel;
  BHmodel->SetLowEnergyUsageLimit(minEn);
  BHmodel->SetHighEnergyUsageLimit(maxEn);
  BHmodel->SetUseSamplingTables(useAlias);
  BHmodel->GetListActiveRegions().resize(1);
  BHmodel->GetListActiveRegions()[0]        = true;
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  BHmodel->Initialize();
  return BHmodel;
}

VecBetheHeitlerPairModel *PrepareVecBHModel(bool useAlias)
{
  VecBetheHeitlerPairModel *BHModel = new VecBetheHeitlerPairModel;
  BHModel->SetLowEnergyUsageLimit(minEn);
  BHModel->SetHighEnergyUsageLimit(maxEn);
  BHModel->SetUseSamplingTables(useAlias);
  BHModel->GetListActiveRegions().resize(1);
  BHModel->GetListActiveRegions()[0] = true;
  BHModel->Initialize();
  return BHModel;
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
    double eKin = minPrimaryEn + (maxEn - minPrimaryEn) * rng.uniform();
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
    double eKin = minPrimaryEn + (maxEn - minPrimaryEn) * rng.uniform();
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

#endif // GEANTV_KLEINNISHINATESTCOMMON_H
