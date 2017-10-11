#include "GeantRunManager.h"

#include "base/Stopwatch.h"
#include "GeantConfig.h"
#include "Geant/Error.h"
#include "VBconnector.h"
#include "GeantPropagator.h"
#include "TaskBroker.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "StdApplication.h"
#include "GeantVTaskMgr.h"
#include "GeantVDetectorConstruction.h"
#include "MCTruthMgr.h"
#include "PrimaryGenerator.h"
#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "LocalityManager.h"

#ifdef USE_ROOT
#include "TApplication.h"
#include "TCanvas.h"
#endif

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#include "Material.h"
#include "Element.h"
#ifdef USE_ROOT
#include "management/RootGeoManager.h"
#endif
#include "volumes/PlacedVolume.h"
#else
#ifdef USE_ROOT
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#endif
#endif

// The classes for integrating in a non-uniform magnetic field
#include "TUniformMagField.h"
#include "FieldEquationFactory.h"
#include "StepperFactory.h"
#include "GUIntegrationDriver.h"

#include "GUFieldPropagator.h"
#include "GUFieldPropagatorPool.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace vecgeom;

//______________________________________________________________________________
GeantRunManager::GeantRunManager(unsigned int npropagators, unsigned int nthreads,
                                 GeantConfig *config)
  : fInitialized(false), fNpropagators(npropagators), fNthreads(nthreads),
    fConfig(config) {
  fPriorityEvents.store(0);
  fTaskId.store(0);
}

//______________________________________________________________________________
bool GeantRunManager::Initialize() {
  // Initialization of run manager
  if (!fNthreads) {
    // Autodiscovery mode using NUMA detection
    fNthreads = 1;   // disabled detection for now
  }

  if (!fNpropagators) {
    Print("Initialize", "Number of propagators set to 1");
    fNpropagators = 1;
  }

  // Not more propagators than events
  if (fNpropagators > fConfig->fNtotal && !fConfig->fUseV3) {
    Print("Initialize", "Number of propagators set to %d", fConfig->fNtotal);
    fNpropagators = fConfig->fNtotal;
  }

  // Increase buffer to give a fair share to each propagator
  int nbuffmax = fConfig->fNtotal/fNpropagators;
  if (fConfig->fUseV3 && nbuffmax == 0)
    nbuffmax = 1;
  if (fConfig->fNbuff > nbuffmax) {
    Print("Initialize", "Number of buffered events reduced to %d", nbuffmax);
    fConfig->fNbuff = nbuffmax;
  }
  fNbuff = fConfig->fNbuff;
  //fConfig->fNbuff *= fNpropagators;
  fConfig->fMaxPerEvent = 5 * fConfig->fNaverage;
  fConfig->fMaxTracks = fConfig->fMaxPerEvent * fConfig->fNbuff;

  if (!fPrimaryGenerator) {
    Fatal("GeantRunManager::Initialize", "The primary generator has to be defined");
    return false;
  }

  if (!fApplication) {
    Fatal("GeantRunManager::Initialize", "The user application has to be defined");
    return false;
  }

  if (!fDetConstruction) {
    Warning("GeantRunManager::Initialize", "The user detector construction has to be defined");
    // return false;
  }

//  fPrimaryGenerator->InitPrimaryGenerator();

  for (auto i=0; i<fNpropagators; ++i) {
    GeantPropagator *prop = new GeantPropagator(fNthreads);
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

  // Temporary workaround to allow migration to detector construction
  if (fDetConstruction) {
    fDetConstruction->CreateMaterials();
    fDetConstruction->CreateGeometry();
    fNvolumes = fDetConstruction->SetupGeometry(fVolumes);
  } else {
    LoadGeometry(fConfig->fGeomFileName.c_str());
  }
#ifdef USE_VECGEOM_NAVIGATOR
    fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
#else
    fConfig->fMaxDepth = TGeoManager::GetMaxLevels();
#endif
  Info("GeantRunManager::Initialize", "Geometry created with maxdepth %d\n", fConfig->fMaxDepth);

  // Now we know the geometry depth: create the track data manager
  TrackDataMgr *dataMgr = TrackDataMgr::GetInstance(fConfig->fMaxDepth);

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

  // Configure the locality manager and create the tracks
  LocalityManager *mgr = LocalityManager::Instance();
  if (!mgr->IsInitialized()) {
    mgr->SetNblocks(100);     // <- must be configurable
    mgr->SetBlockSize(1000);  // <- must be configurable
    mgr->Init();
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    if (fConfig->fUseNuma) {
      int nnodes = mgr->GetPolicy().GetNnumaNodes();
      mgr->SetPolicy(NumaPolicy::kCompact);
      // Loop propagatpors and assign NUMA nodes
      if (fNpropagators > 1) {
        for (auto i=0; i<fNpropagators; ++i)
          fPropagators[i]->SetNuma(i % nnodes);
      }
    } else {
      mgr->SetPolicy(NumaPolicy::kSysDefault);
    }
#endif
  }

  fDoneEvents = BitSet::MakeInstance(fConfig->fNtotal);

  if (fConfig->fUseRungeKutta) {
    PrepareRkIntegration();
  }

  fPrimaryGenerator->InitPrimaryGenerator();

  // Initialize task data
  int nthreads = GetNthreadsTotal();
  fTDManager = new TDManager(nthreads, fConfig->fMaxPerBasket);

  // Initialize application
  if (fConfig->fUseStdScoring) {
    fStdApplication = new StdApplication(this);
    fStdApplication->Initialize();
  }
  fApplication->Initialize();

  // Attach user data and physics data to task data
  for (int i = 0; i < nthreads; i++) {
    GeantTaskData *td = fTDManager->GetTaskData(i);
    if (fPhysicsInterface) fPhysicsInterface->AttachUserData(td);
    if (fStdApplication) fStdApplication->AttachUserData(td);
    fApplication->AttachUserData(td);
  }

  // Initialize all propagators
  for (auto i=0; i<fNpropagators; ++i)
    fPropagators[i]->Initialize();

  GeantTaskData *td = fTDManager->GetTaskData(0);
  td->AttachPropagator(fPropagators[0], 0);

  // Initialize the event server
  fEventServer = new GeantEventServer(fConfig->fNbuff, this);
  
  dataMgr->Print();
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
  delete fPhysicsInterface;
  delete fVectorPhysicsProcess;
  if (fStdApplication) {
    fStdApplication->DeleteUserData();
    delete fStdApplication;
  }
  fApplication->DeleteUserData();
  delete fApplication;
  delete fTaskMgr;

  delete fTDManager;
//  if (fTaskData) {
//    for (auto i = 0; i < fNthreads; i++)
//      delete fTaskData[i];
//    delete[] fTaskData;
//  }

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
    GeantVDetectorConstruction::LoadVecGeomGeometry(fBroker);
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
#else
    geom->SetMaxThreads(GetNthreadsTotal());
    TObjArray *lvolumes = geom->GetListOfVolumes();
    fNvolumes = lvolumes->GetEntries();
    for (auto ivol = 0; ivol < fNvolumes; ivol++)
      fVolumes.push_back((TGeoVolume *)lvolumes->At(ivol));
#endif
  } else {
    Error("GeantPropagator::LoadGeometry", "Cannot load geometry from file %s", filename);
    return false;
  }
#else
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    geom->LoadGeometryFromSharedLib(filename);
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
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
void GeantRunManager::EventTransported(GeantEvent *event, GeantTaskData *td)
{
// Actions executed after an event is transported.
  // Adjust number of prioritized events
  if (event->IsPrioritized()) fPriorityEvents--;
  // closing event in MCTruthManager
  if(fTruthMgr) fTruthMgr->CloseEvent(event->GetEvent());
  // event->Print();
  // Digitizer
  Info("EventTransported", " = task %d digitizing event %d with %d tracks", td->fTid, event->GetEvent(), event->GetNtracks());
//  LocalityManager *lmgr = LocalityManager::Instance();
//  Printf("   NQUEUED = %d  NBLOCKS = %d NRELEASED = %d",
//         lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
  fApplication->FinishEvent(event);
  // Signal completion of one event to the event server
  fDoneEvents->SetBitNumber(event->GetEvent());
  assert(event->GetNtracks() > 0);
  fEventServer->CompletedEvent(event, td);
}

//______________________________________________________________________________
int GeantRunManager::ProvideWorkTo(GeantPropagator *prop)
{
// Provide work to a given propagator which became idle
  int nshared = 0;
  for (auto i=0; i<fNpropagators; ++i) {
    if (fPropagators[i] == prop) continue;
    nshared += fPropagators[i]->ShareWork(*prop);
  }
  // if (nshared) Printf("Propagator %p stole %d baskets", prop, nshared);
  return nshared;
}

//______________________________________________________________________________
GeantPropagator *GeantRunManager::GetIdlePropagator() const {
// Returns the first found idle propagator if any
  for (auto i=0; i<fNpropagators; ++i) {
    if (fPropagators[i]->fCompleted) continue;
    if (fPropagators[i]->IsIdle()) return fPropagators[i];
  }
  return nullptr;
}

//______________________________________________________________________________
void GeantRunManager::RunSimulation() {
  // Start simulation for all propagators
  Initialize();

  Printf("==========================================================================");
  Printf("= GeantV run started with %d propagator(s) using %d worker threads each ====", fNpropagators, fNthreads);
  if (fConfig->fUsePhysics)
    Printf("  Physics ON");
  else
    Printf("  Physics OFF");
  if (fConfig->fUseRungeKutta)
    Printf("  Runge-Kutta integration ON with epsilon= %g", fConfig->fEpsilonRK);
  else
    Printf("  Runge-Kutta integration OFF");
  Printf("==========================================================================");
#ifdef USE_ROOT
  if (fConfig->fUseMonitoring)
    new TCanvas("cscheduler", "Scheduler monitor", 900, 600);
  if (fConfig->fUseAppMonitoring)
    new TCanvas("capp", "Application canvas", 700, 800);
#endif
  vecgeom::Stopwatch timer;
  timer.Start();
  for (auto i=0; i<fNpropagators; ++i)
    fListThreads.emplace_back(GeantPropagator::RunSimulation, fPropagators[i], fNthreads);

  for (auto &t : fListThreads) {
    t.join();
  }
  timer.Stop();
  double rtime = timer.Elapsed();
  double ctime = timer.CpuElapsed();
  long ntransported = fEventServer->GetNprimaries();
  long nsteps = 0;
  long nsnext = 0;
  long nphys = 0;
  long nmag = 0;
  long nsmall = 0;
  long ncross = 0;
  long npushed = 0;
  long nkilled = 0;

  for (auto i=0; i<fNpropagators; ++i) {
    ntransported += fPropagators[i]->fNtransported.load();
    nsteps += fPropagators[i]->fNsteps.load();
    nsnext += fPropagators[i]->fNsnext.load();
    nphys += fPropagators[i]->fNphys.load();
    nmag += fPropagators[i]->fNmag.load();
    nsmall += fPropagators[i]->fNsmall.load();
    ncross += fPropagators[i]->fNcross.load();
    npushed += fPropagators[i]->fNpushed.load();
    nkilled += fPropagators[i]->fNkilled.load();
  }
  Printf("=== Summary: %d propagators x %d threads: %ld primaries/%ld tracks,  total steps: %ld, snext calls: %ld, "
         "phys steps: %ld, mag. field steps: %ld, small steps: %ld, pushed: %ld, killed: %ld, bdr. crossings: %ld  RealTime=%gs CpuTime=%gs",
         fNpropagators, fNthreads, fEventServer->GetNprimaries(), ntransported, nsteps, nsnext, nphys, nmag, nsmall, npushed, nkilled, ncross, rtime, ctime);
  //LocalityManager *lmgr = LocalityManager::Instance();
  // Printf("NQUEUED = %d  NBLOCKS = %d NRELEASED = %d",
  //       lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
#ifdef USE_VECGEOM_NAVIGATOR
  Printf("=== Navigation done using VecGeom ====");
#else
  Printf("=== Navigation done using TGeo    ====");
#endif

  FinishRun();
#ifdef USE_ROOT
  if (gApplication) delete gApplication;
#endif
}

//______________________________________________________________________________
bool GeantRunManager::FinishRun() {
  // Run termination actions.
  if (fTaskMgr) fTaskMgr->Finalize();
  fApplication->FinishRun();
  if (fStdApplication)
    fStdApplication->FinishRun();
  // Actions to follow
  return true;
}

//______________________________________________________________________________
void GeantRunManager::StopTransport() {
  // Signal all propagators that transport has stopped
  for (auto i=0; i<fNpropagators; ++i) {
    fPropagators[i]->StopTransport();
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
