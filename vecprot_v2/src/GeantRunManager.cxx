#include "GeantRunManager.h"

#include "GeantConfig.h"
#include "Geant/Error.h"
#include "GeantPropagator.h"
#include "TaskBroker.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#ifdef USE_ROOT
#include "management/RootGeoManager.h"
#endif
#include "volumes/PlacedVolume.h"
#else
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#endif

using namespace Geant;
using namespace vecgeom;

//______________________________________________________________________________
GeantRunManager::GeantRunManager(unsigned int npropagators, unsigned int nthreads,
                                 GeantConfig *config) 
  : fInitialized(false), fNpropagators(npropagators), fNthreads(nthreads),
    fConfig(config) {
}

//______________________________________________________________________________
bool GeantRunManager::Initialize() {
  // Initialization of run manager
  if (!fNthreads) {
    // Autodiscovery mode using NUMA detection
    fNthreads = 1;   // disabled detection for now
  }
  
  if (!fNpropagators) {
    fNpropagators = 1;
  }
  for (auto i=0; i<fNpropagators; ++i) {
    GeantPropagator *prop = GeantPropagator::NewInstance(fNthreads);
    fPropagators.push_back(prop);
    prop->fRunMgr = this;
    prop->SetConfig(fConfig);
    prop->fApplication = fApplication;
    prop->fStdApplication = fStdApplication;
    prop->fTaskMgr = fTaskMgr;
    prop->fProcess = fProcess;
    prop->fVectorPhysicsProcess = fVectorPhysicsProcess;
    prop->fPrimaryGenerator = fPrimaryGenerator;
    prop->fTruthMgr = fTruthMgr;
  }
  
#ifndef VECCORE_CUDA
  LoadGeometry(fConfig->fGeomFileName.c_str());
#endif
  fInitialized = true;
  return fInitialized;
}

//______________________________________________________________________________
GeantRunManager::~GeantRunManager() {
  for (auto i=0; i<fNpropagators; ++i) delete fPropagators[i];
  fPropagators.clear();
  for (auto i=0; i<fNvolumes; ++i) {
    Volume_t *vol = (Volume_t*)fVolumes[i];
#ifdef USE_VECGEOM_NAVIGATOR
    VBconnector *connector = (VBconnector*)vol->GetBasketManagerPtr();
    vol->SetBasketManagerPtr(nullptr);
#else
    VBconnector *connector = (VBconnector*)vol->GetFWExtension();
    vol->SetFWExtension(nullptr);
#endif
    delete connector;
  }
  delete fConfig;
}

//______________________________________________________________________________
bool GeantRunManager::LoadGeometry(const char *filename) {
// Load geometry from given file.
#ifdef USE_ROOT
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
#else
  TGeoManager *geom = (gGeoManager) ? gGeoManager : TGeoManager::Import(filename);
#endif
  if (geom) {
#ifdef USE_VECGEOM_NAVIGATOR
    LoadVecGeomGeometry();
    fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
#else
    fConfig->fMaxDepth = TGeoManager::GetMaxLevels();
    TObjArray *lvolumes = gGeoManager->GetListOfVolumes();
    fNvolumes = lvolumes->GetEntries();
    for (auto ivol = 0; ivol < fNvolumes; ivol++)
      fVolumes.push_back((TGeoVolume *)lvolumes->At(ivol));
#endif
    printf("Geometry depth %d\n", fConfig->fMaxDepth);
    return true;
  }
  Geant::Error("GeantPropagator::LoadGeometry", "Cannot load geometry from file %s", filename);
  return false;
#else
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    geom->LoadGeometryFromSharedLib(filename);
    fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
    return true;
  }
  return false;
#endif
  for (auto i=0; i<fNvolumes; ++i) {
    Volume_t *vol = (Volume_t*)fVolumes[i];
    VBconnector *connector = new VBconnector(i);
#ifdef USE_VECGEOM_NAVIGATOR
    vol->SetBasketManagerPtr(connector);
#else
    vol->SetFWExtension(connector);
#endif
  } 
}  

//______________________________________________________________________________
bool GeantRunManager::LoadVecGeomGeometry() {
#ifdef USE_VECGEOM_NAVIGATOR
  if (vecgeom::GeoManager::Instance().GetWorld() == NULL) {
    printf("Now loading VecGeom geometry\n");
    vecgeom::RootGeoManager::Instance().LoadRootGeometry();
    printf("Loading VecGeom geometry done\n");
    std::vector<vecgeom::LogicalVolume *> v1;
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(v1);
    printf("Have logical volumes %ld\n", v1.size());
    std::vector<vecgeom::VPlacedVolume *> v2;
    vecgeom::GeoManager::Instance().getAllPlacedVolumes(v2);
    printf("Have placed volumes %ld\n", v2.size());
    //    vecgeom::RootGeoManager::Instance().world()->PrintContent();

  }
  if (fBroker) {
    printf("Now upload VecGeom geometry to Coprocessor(s)\n");
    return fBroker->UploadGeometry();
  }
  InitNavigators();
  return true;
#else
  return false;
#endif
}

//______________________________________________________________________________
void GeantRunManager::InitNavigators() {
#if USE_VECGEOM_NAVIGATOR == 1
  for (auto &lvol : GeoManager::Instance().GetLogicalVolumesMap()) {
    if (lvol.second->GetDaughtersp()->size() < 4) {
      lvol.second->SetNavigator(NewSimpleNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 5) {
      lvol.second->SetNavigator(SimpleABBoxNavigator<>::Instance());
    }
    if (lvol.second->GetDaughtersp()->size() >= 10) {
      lvol.second->SetNavigator(HybridNavigator<>::Instance());
      HybridManager2::Instance().InitStructure((lvol.second));
    }
    lvol.second->SetLevelLocator(SimpleABBoxLevelLocator::GetInstance());
  }
#endif
}

//______________________________________________________________________________
void GeantRunManager::RunSimulation() {
  // Start simulation for all propagators
  Initialize();
  // Should start threads instead
  for (auto i=0; i<fNpropagators; ++i)
    fPropagators[i]->PropagatorGeom("", fNthreads);
  FinishRun();
}

//______________________________________________________________________________
bool GeantRunManager::FinishRun() {
  // Run termination actions.
  for (auto i=0; i<fNpropagators; ++i) {
    fNprimaries += fPropagators[i]->GetNprimaries();
  }
  // Actions to follow
  return true;
}
