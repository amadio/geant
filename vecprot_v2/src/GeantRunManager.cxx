#include "GeantRunManager.h"

#include "base/Stopwatch.h"
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
#include "GeantEvent.h"
#include "GeantEventServer.h"
#include "LocalityManager.h"

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
  if (fNpropagators > fConfig->fNtotal) {
    Print("Initialize", "Number of propagators set to %d", fConfig->fNtotal);
    fNpropagators = fConfig->fNtotal;
  }

  // Increase buffer to give a fair share to each propagator
  int nbuffmax = fConfig->fNtotal/fNpropagators;
  if (fConfig->fNbuff > nbuffmax) {
    Print("Initialize", "Number of buffered events reduced to %d", nbuffmax);
    fConfig->fNbuff = nbuffmax;
  }
  fNbuff = fConfig->fNbuff;
  fConfig->fNbuff *= fNpropagators;
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
  fPrimaryGenerator->InitPrimaryGenerator();
  fEventServer = new GeantEventServer(fConfig->fNtotal, this);
  for (int i=0; i<fConfig->fNtotal; ++i)
    fEventServer->AddEvent();

  for (auto i=0; i<fNpropagators; ++i)
    fPropagators[i]->Initialize();

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
  delete fApplication;
  delete fTaskMgr;

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
    geom->SetMaxThreads(GetNthreadsTotal());
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
#ifdef USE_ROOT
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
#endif
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
void GeantRunManager::EventTransported(int evt)
{
// Actions executed after an event is transported.
  // Signal completion of one event to the event server
  fEventServer->CompletedEvent(evt);
  // Adjust number of prioritized events
  GeantEvent *event = fEventServer->GetEvent(evt);
  if (event->IsPrioritized()) fPriorityEvents--;
  // closing event in MCTruthManager
  if(fTruthMgr) fTruthMgr->CloseEvent(evt);
  event->Print();
  // Digitizer (todo)
  Info("EventTransported", " = digitizing event %d with %d tracks", evt, event->GetNtracks());
  LocalityManager *lmgr = LocalityManager::Instance();
  Printf("   NQUEUED = %d  NBLOCKS = %d NRELEASED = %d", 
         lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
  //fApplication->Digitize(event);
  fDoneEvents->SetBitNumber(evt);
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
  long ntransported = fNprimaries;
  long nsteps = 0;
  long nsnext = 0;
  long nphys = 0;
  long nmag = 0;
  long nsmall = 0;
  long ncross = 0;

  for (auto i=0; i<fNpropagators; ++i) {
    ntransported += fPropagators[i]->fNtransported.load();
    nsteps += fPropagators[i]->fNsteps.load();
    nsnext += fPropagators[i]->fNsnext.load();
    nphys += fPropagators[i]->fNphys.load();
    nmag += fPropagators[i]->fNmag.load();
    nsmall += fPropagators[i]->fNsmall.load();
    ncross += fPropagators[i]->fNcross.load();
  }
  Printf("=== Summary: %d propagators x %d threads: %ld primaries/%ld tracks,  total steps: %ld, snext calls: %ld, "
         "phys steps: %ld, mag. field steps: %ld, small steps: %ld bdr. crossings: %ld  RealTime=%gs CpuTime=%gs",
         fNpropagators, fNthreads, fNprimaries, ntransported, nsteps, nsnext, nphys, nmag, nsmall, ncross, rtime, ctime);
  LocalityManager *lmgr = LocalityManager::Instance();
  Printf("NQUEUED = %d  NBLOCKS = %d NRELEASED = %d", 
         lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
#ifdef USE_VECGEOM_NAVIGATOR
  Printf("=== Navigation done using VecGeom ====");
#else
  Printf("=== Navigation done using TGeo    ====");
#endif

  FinishRun();
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
