#include "GeantRunManager.h"

#include "GeantConfig.h"
#include "Geant/Error.h"
#include "GeantPropagator.h"
#include "TaskBroker.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "StdApplication.h"
#include "GeantVTaskMgr.h"
#include "MCTruthMgr.h"
#include "PrimaryGenerator.h"

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

// The classes for integrating in a non-uniform magnetic field
#include "TUniformMagField.h"
#include "FieldEquationFactory.h"
#include "StepperFactory.h"
#include "GUIntegrationDriver.h"

#include "GUFieldPropagator.h"
#include "GUFieldPropagatorPool.h"

using namespace Geant;
using namespace vecgeom;

//______________________________________________________________________________
GeantRunManager::GeantRunManager(unsigned int npropagators, unsigned int nthreads,
                                 GeantConfig *config) 
  : fInitialized(false), fNpropagators(npropagators), fNthreads(nthreads),
    fConfig(config) {
  fFeederLock.clear();
  fPriorityEvents.store(0);
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

  fConfig->fMaxPerEvent = 5 * fConfig->fNaverage;
  fConfig->fMaxTracks = fConfig->fMaxPerEvent * fConfig->fNbuff;
  
  for (auto i=0; i<fNpropagators; ++i) {
    GeantPropagator *prop = GeantPropagator::NewInstance(fNthreads);
    fPropagators.push_back(prop);
    prop->fRunMgr = this;
    prop->SetConfig(fConfig);
    prop->fApplication = fApplication;
    prop->fStdApplication = fStdApplication;
    prop->fTaskMgr = fTaskMgr;
    prop->fProcess = fProcess;
    prop->fPhysicsInterface = fPhysicsInterface;
    prop->fVectorPhysicsProcess = fVectorPhysicsProcess;
    prop->fPrimaryGenerator = fPrimaryGenerator;
    prop->fTruthMgr = fTruthMgr;
  }
  
//#ifndef VECCORE_CUDA
  LoadGeometry(fConfig->fGeomFileName.c_str());
//#endif

  fDoneEvents = BitSet::MakeInstance(fConfig->fNtotal);

  // Initialize the process(es)
#ifdef USE_REAL_PHYSICS
  if (!fPhysicsInterface) {
    Geant::Fatal("GeantRunManager::Initialize", "The physics process interface has to be initialized before this");
    return false;
  }
  // Initialize the physics
  fPhysicsInterface->Initialize();
#else
  if (!fProcess) {
    Geant::Fatal("GeantRunManager::Initialize", "The physics process has to be initialized before this");
    return false;
  }
  // Initialize the process(es)
  fProcess->Initialize();
  #if USE_VECPHYS == 1
  fVectorPhysicsProcess->Initialize();
  #endif
#endif

  if (fConfig->fUseRungeKutta) {
    PrepareRkIntegration();
  }

  fNtracks = new int[fConfig->fNbuff];
  memset(fNtracks, 0, fConfig->fNbuff * sizeof(int));

  int nthreads = GetNthreadsTotal();
  fTaskData = new GeantTaskData *[nthreads];
  for (int i = 0; i < nthreads; i++) {
    fTaskData[i] = new GeantTaskData(nthreads, fConfig->fMaxDepth, fConfig->fMaxPerBasket);
    fTaskData[i]->fTid = i;
  }
  if (fConfig->fUseStdScoring) {
    fStdApplication = new StdApplication(this);
    fStdApplication->Initialize();
  }
  fApplication->Initialize();

  for (auto i=0; i<fNpropagators; ++i)
    fPropagators[i]->Initialize();

  fEvents = new GeantEvent *[fConfig->fNbuff];
  memset(fEvents, 0, fConfig->fNbuff * sizeof(GeantEvent *));

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
  BitSet::ReleaseInstance(fDoneEvents);
  delete fProcess;
  delete fVectorPhysicsProcess;
  delete fApplication;
  delete fTaskMgr;
  if (fEvents) {
    for (auto i = 0; i < fConfig->fNbuff; i++)
      delete fEvents[i];
    delete[] fEvents;
  }

  if (fTaskData) {
    for (auto i = 0; i < fNthreads; i++)
      delete fTaskData[i];
    delete[] fTaskData;
  }

  delete fConfig;
}

//______________________________________________________________________________
bool GeantRunManager::LoadGeometry(const char *filename) {
// Load geometry from given file.
#ifdef USE_ROOT
  if (!gGeoManager) TGeoManager::Import(filename);
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
#else
  TGeoManager *geom = gGeoManager;
#endif
  if (geom) {
#ifdef USE_VECGEOM_NAVIGATOR
    LoadVecGeomGeometry();
    fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
#else
    fConfig->fMaxDepth = TGeoManager::GetMaxLevels();
    TObjArray *lvolumes = geom->GetListOfVolumes();
    fNvolumes = lvolumes->GetEntries();
    for (auto ivol = 0; ivol < fNvolumes; ivol++)
      fVolumes.push_back((TGeoVolume *)lvolumes->At(ivol));
#endif
    Info("GeantRunManager::LoadGeometry", "Geometry depth %d\n", fConfig->fMaxDepth);
  } else {
    Error("GeantPropagator::LoadGeometry", "Cannot load geometry from file %s", filename);
    return false;
  }
#else
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    geom->LoadGeometryFromSharedLib(filename);
    fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
    Info("GeantRunManager::LoadGeometry", "Geometry depth %d\n", fConfig->fMaxDepth);
  } else {
    Error("GeantPropagator::LoadGeometry", "Cannot load geometry from file %s", filename);
    return false;
  }
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
  return true;
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
void GeantRunManager::PrepareRkIntegration() {

  using GUFieldPropagatorPool = ::GUFieldPropagatorPool;
  using GUFieldPropagator = ::GUFieldPropagator;

  // Initialise the classes required for tracking in field
  const unsigned int Nvar = 6; // Integration will occur over 3-position & 3-momentum coord.
  using Field_t = TUniformMagField;
  using Equation_t = TMagFieldEquation<Field_t, Nvar>;

  auto gvField = new Field_t(fieldUnits::kilogauss * ThreeVector(0.0, 0.0, fConfig->fBmag));
  auto gvEquation = FieldEquationFactory::CreateMagEquation<Field_t>(gvField);

  GUVIntegrationStepper *aStepper = StepperFactory::CreateStepper<Equation_t>(gvEquation); // Default stepper

  const double hminimum = 1.0e-5; // * centimeter; =  0.0001 * millimeter;  // Minimum step = 0.1 microns
  // const double epsTol = 3.0e-4;               // Relative error tolerance of integration
  int statisticsVerbosity = 0;
  cout << "Parameters for RK integration in magnetic field: " << endl;
  cout << "   Driver parameters:  eps_tol= " << fConfig->fEpsilonRK << "  h_min= " << hminimum << endl;

  auto integrDriver = new GUIntegrationDriver(hminimum, aStepper, Nvar, statisticsVerbosity);
  // GUFieldPropagator *
  auto fieldPropagator = new GUFieldPropagator(integrDriver, fConfig->fEpsilonRK); // epsTol);

  static GUFieldPropagatorPool *fpPool = GUFieldPropagatorPool::Instance();
  assert(fpPool); // Cannot be zero
  if (fpPool) {
    fpPool->RegisterPrototype(fieldPropagator);
    // Create clones for other threads
    fpPool->Initialize(fNthreads);
  } else {
    Geant::Error("PrepareRkIntegration", "Cannot find GUFieldPropagatorPool Instance.");
  }
}

//______________________________________________________________________________
int GeantRunManager::Feeder(GeantTaskData *td) {
  // Feeder called by any thread to inject the next event(s)
  // Only one thread at a time
  if (fFeederLock.test_and_set(std::memory_order_acquire))
    return -1;
  int nbaskets = 0;
  if (!fConfig->fLastEvent) {
    nbaskets = ImportTracks(fConfig->fNbuff, 0, 0, td);
    fConfig->fLastEvent = fConfig->fNbuff;
    fFeederLock.clear(std::memory_order_release);
    return nbaskets;
  }
  // Check and mark finished events
  for (int islot = 0; islot < fConfig->fNbuff; islot++) {
    GeantEvent *evt = fEvents[islot];
    if (fDoneEvents->TestBitNumber(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      fPriorityEvents--;
      evt->Print();
      
      // closing event in MCTruthManager
      if(fTruthMgr) fTruthMgr->CloseEvent(evt->GetEvent());
      
      // Digitizer (todo)
      int ntracks = fNtracks[islot];
      Info("Feeder:", " = digitizing event %d with %d tracks", evt->GetEvent(), ntracks);
      //            propagator->fApplication->Digitize(evt->GetEvent());
      fDoneEvents->SetBitNumber(evt->GetEvent());
      if (fConfig->fLastEvent < fConfig->fNtotal) {
      Info("Feeder:", "  => Importing event %d", fConfig->fLastEvent);
        nbaskets += ImportTracks(1, fConfig->fLastEvent, islot, td);
        fConfig->fLastEvent++;
      }
    }
  }

  fFeederLock.clear(std::memory_order_release);
  return nbaskets;
}

//______________________________________________________________________________
int GeantRunManager::ImportTracks(int nevents, int startevent, int startslot, GeantTaskData *thread_data) {
// Import tracks from "somewhere". Here we just generate nevents.

#ifdef USE_VECGEOM_NAVIGATOR
  using vecgeom::SimpleNavigator;
  using vecgeom::Vector3D;
  using vecgeom::Precision;
  using vecgeom::GeoManager;
#endif

  GeantPropagator *prop = thread_data->fPropagator;
  Volume_t *vol = 0;
  int ntracks = 0;
  int ntotal = 0;
  int ndispatched = 0;
  GeantTaskData *td = thread_data;
  VolumePath_t *startpath = VolumePath_t::MakeInstance(fConfig->fMaxDepth);

#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::SimpleNavigator nav;
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
#endif

  static bool init = true;
  if (init)
    init = false;
  int event = startevent;
  GeantEventInfo eventinfo;
  for (int slot = startslot; slot < startslot + nevents; slot++) {
    eventinfo = fPrimaryGenerator->NextEvent();
// Set initial track states
#ifdef USE_VECGEOM_NAVIGATOR
    startpath->Clear();
    nav.LocatePoint(GeoManager::Instance().GetWorld(),
                    Vector3D<Precision>(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert), *startpath, true);
    vol = const_cast<Volume_t *>(startpath->Top()->GetLogicalVolume());
#else
    TGeoNode *node = nav->FindNode(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);
    vol = node->GetVolume();
    startpath->InitFromNavigator(nav);
#endif
    td->fVolume = vol;
    ntracks = eventinfo.ntracks;
    ntotal += ntracks;
    fNprimaries += ntracks;
    if (!fEvents[slot])
      fEvents[slot] = new GeantEvent();
    fEvents[slot]->SetSlot(slot);
    fEvents[slot]->SetEvent(event);
    fEvents[slot]->Reset();
    // Set priority threshold to non-default value
    if (fConfig->fPriorityThr > 0)
      fEvents[slot]->SetPriorityThr(fConfig->fPriorityThr);

    // start new event in MCTruthMgr
    if(fTruthMgr) fTruthMgr->OpenEvent(fEvents[slot]->GetEvent());

    for (int i = 0; i < ntracks; i++) {
      GeantTrack &track = td->GetTrack();
      track.SetPath(startpath);
      track.SetNextPath(startpath);
      track.SetEvent(event);
      track.SetEvslot(slot);
      fPrimaryGenerator->GetTrack(i, track);
      if (!track.IsNormalized())
        track.Print("Not normalized");
      track.fBoundary = false;
      track.fStatus = kAlive;
      prop->AddTrack(track);
      ndispatched += prop->DispatchTrack(track, thread_data);
    }
    event++;
  }

  VolumePath_t::ReleaseInstance(startpath);
  int tid = td->fTid;
  Geant::Print("ImportTracks", "[%d] Imported %d tracks from events %d to %d. Dispatched %d baskets.", tid, ntotal, startevent,
               startevent + nevents - 1, ndispatched);
  return ndispatched;
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
