#ifndef GEANTV_SAUTERGAVRILATESTCOMMON_H
#define GEANTV_SAUTERGAVRILATESTCOMMON_H

#include "Geant/LightTrack.h"
#include "Geant/TaskData.h"
#include "Geant/SauterGavrilaPhotoElectricModel.h"
#include "Common.h"
#include <random>
#include "Geant/MaterialCuts.h"
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"
#include "Geant/Region.h"
#include "Geant/PhysicsParameters.h"
#include "Geant/Particle.h"
#include "Geant/Electron.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"

using geantphysics::LightTrack;
using geantphysics::LightTrack_v;
using geantphysics::SauterGavrilaPhotoElectricModel;
using geantphysics::Material;
using geantphysics::MaterialCuts;

const int kMinBasket = 256;
const int kMaxBasket = 256;
const double minEn   = 100 * geant::units::eV;
const double maxEn   = 50 * geant::units::keV;

SauterGavrilaPhotoElectricModel *PrepareSauterGavrilaModel(bool useAlias)
{
  SauterGavrilaPhotoElectricModel *sauterG = new SauterGavrilaPhotoElectricModel();
  sauterG->SetLowEnergyUsageLimit(minEn);
  sauterG->SetHighEnergyUsageLimit(maxEn);
  sauterG->SetUseSamplingTables(useAlias);
  (sauterG->GetListActiveRegions()).resize(1); // one region
  (sauterG->GetListActiveRegions())[0] = true; // make it active there
  // - Create one PhysicsParameters object (with defult values)
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  // - Initialisation of the model
  sauterG->Initialize();
  return sauterG;
}

SauterGavrilaPhotoElectricModel *PrepareVecSauterGavrilaModel(bool useAlias)
{
  SauterGavrilaPhotoElectricModel *sauterG = new SauterGavrilaPhotoElectricModel();
  sauterG->SetLowEnergyUsageLimit(minEn);
  sauterG->SetHighEnergyUsageLimit(maxEn);
  sauterG->SetUseSamplingTables(useAlias);
  (sauterG->GetListActiveRegions()).resize(1); // one region
  (sauterG->GetListActiveRegions())[0] = true; // make it active there
  // - Create one PhysicsParameters object (with defult values)
  geantphysics::PhysicsParameters *physPars = new geantphysics::PhysicsParameters();
  // - Set it to be active in region index 0
  (physPars->GetListActiveRegions()).resize(1);
  (physPars->GetListActiveRegions())[0] = true;
  sauterG->Initialize();
  return sauterG;
}

void PreparePrimaries(std::vector<LightTrack> &output, int N)
{
  geant::RngWrapper rng;
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
    gamma.SetMaterialCutCoupleIndex(0);
    output.push_back(gamma);
  }
}

void PreparePrimaries(LightTrack_v &output, int N)
{
  geant::RngWrapper rng;
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
    output.SetMaterialCutCoupleIndex(0,i);
  }
}

void SetUpSimulation()
{
    vecgeom::GeoManager::Instance().Clear();
    static std::string   particleName("gamma");
    static std::string   materialName("NIST_MAT_CONCRETE"); // material is CONCRETE
    static std::string   photoElectricModelName("SauterGavrilaPhotoElectric");
    //static double        prodCutValue      = 0.1;             // by default in length and internal units i.e. [cm]
    //static bool          isProdCutInLength = true;            // is the production cut value given in length ?
    
    //============================= Set user defined input data =================================//
    // Create target material: which is supposed to be a NIST Material
    Material *matDetector = Material::NISTMaterial(materialName);
    
    // Set production cuts if needed: not needed in case of pair-production
    bool iscutinlength  = true;
    double prodCutValue = 1. * geant::units::mm;
    double gcut         = prodCutValue;
    double emcut        = prodCutValue;
    double epcut        = prodCutValue;
    
    // Create primary particle
    //Particle *
    //particle = Gamma::Definition();
    //std::string  pname;
    
    //============= Initialization i.e. building up and init the physics ========================//
    // Create a dummy vecgeom::geometry:
    //  - with only one volume i.e. world
    //  - create a region and set production cuts
    //
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
    
    // if primary particle energy < gamma production cut => there is no secondary gamma production
    // So get the MaterialCuts of the target: we have only one
    //const MaterialCuts *
    //const MaterialCuts *matCut = MaterialCuts::GetMaterialCut(aRegion->GetIndex(),matDetector->GetIndex());
    //===========================================================================================//
}

#endif // GEANTV_SAUTERGAVRILATESTCOMMON_H
