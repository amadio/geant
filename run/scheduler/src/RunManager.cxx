#include "Geant/RunManager.h"

#include "base/Stopwatch.h"
#include "GeantConfig.h"
#include "Geant/Error.h"
#include "Geant/VBconnector.h"
#include "Geant/Propagator.h"
#include "Geant/WorkloadManager.h"
#include "Geant/TaskBroker.h"
#include "Geant/PhysicsInterface.h"
#include "Geant/StdApplication.h"
#include "Geant/UserDetectorConstruction.h"
#include "Geant/MCTruthMgr.h"
#include "Geant/PrimaryGenerator.h"
#include "Geant/Event.h"
#include "Geant/EventSet.h"
#include "Geant/EventServer.h"
#include "Geant/LocalityManager.h"
#include "Geant/SimulationStage.h"
#include "Geant/BasketCounters.h"

#ifdef USE_ROOT
#include "TApplication.h"
#include "TCanvas.h"
#include "management/RootGeoManager.h"
#endif

#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/NewSimpleNavigator.h"
#include "navigation/SimpleABBoxNavigator.h"
#include "navigation/SimpleABBoxLevelLocator.h"
#include "navigation/HybridNavigator2.h"
#include "Geant/Material.h"
#include "Geant/Element.h"
#include "volumes/PlacedVolume.h"

#include "Geant/UserFieldConstruction.h"

// The classes for integrating in a non-uniform magnetic field
#include "Geant/ScalarUniformMagField.h"
#include "Geant/FieldEquationFactory.h"
#include "Geant/StepperFactory.h"
#include "Geant/ScalarIntegrationDriver.h"

#include "Geant/GUFieldPropagator.h"
#include "Geant/GUFieldPropagatorPool.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace vecgeom;

//______________________________________________________________________________
RunManager::RunManager(unsigned int npropagators, unsigned int nthreads, GeantConfig *config)
    : fInitialized(false), fNpropagators(npropagators), fNthreads(nthreads), fConfig(config)
{
  fPriorityEvents.store(0);
  fTaskId.store(0);
  fEventSetsLock.clear();
}

//______________________________________________________________________________
bool RunManager::Initialize()
{
  const char *methodName = "RunManager::Initialize";
  // Initialization of run manager
  if (fInitialized) return fInitialized;
  if (!fNthreads) {
    // Autodiscovery mode using NUMA detection
    fNthreads = 1; // disabled detection for now
  }

  if (!fNpropagators) {
    Print(methodName, "Number of propagators set to 1");
    fNpropagators = 1;
  }

  // Not more propagators than events
  if (fNpropagators > fConfig->fNtotal && !fConfig->fUseV3) {
    Print("Initialize", "Number of propagators set to %d", fConfig->fNtotal);
    fNpropagators = fConfig->fNtotal;
  }

  // Increase buffer to give a fair share to each propagator
  int nbuffmax                                   = fConfig->fNtotal / fNpropagators;
  if (fConfig->fUseV3 && nbuffmax == 0) nbuffmax = 1;
  if (fConfig->fNbuff > nbuffmax) {
    Print("Initialize", "Number of buffered events reduced to %d", nbuffmax);
    fConfig->fNbuff = nbuffmax;
  }
  fNbuff = fConfig->fNbuff;
  // fConfig->fNbuff *= fNpropagators;
  fConfig->fMaxPerEvent = 5 * fConfig->fNaverage;
  fConfig->fMaxTracks   = fConfig->fMaxPerEvent * fConfig->fNbuff;

  bool externalLoop = fConfig->fRunMode == GeantConfig::kExternalLoop;

  if (!fPrimaryGenerator && !externalLoop) {
    Fatal(methodName, "The primary generator has to be defined");
    return false;
  }

  if (!fApplication) {
    Fatal(methodName, "The user application has to be defined");
    return false;
  }

  if (!fDetConstruction) {
    Warning(methodName, "The user detector construction should be defined");
    // return false;
  } else {
    if (!(fDetConstruction->GetFieldConstruction())) {
      Printf("- RunManager::Initialize - %s.\n", " no User Field Construction found in Detector Construction.");
      Printf("    Created a default field= %f %f %f\n", fBfieldArr[0], fBfieldArr[1], fBfieldArr[2]);
      auto fieldCtion = new UserFieldConstruction();
      fieldCtion->UseConstantMagField(fBfieldArr);
      SetUserFieldConstruction(fieldCtion);
      Warning(methodName, "A default field construction was defined");
    }
  }

  //  fPrimaryGenerator->InitPrimaryGenerator();

  for (auto i = 0; i < fNpropagators; ++i) {
    Propagator *prop = new Propagator(fNthreads);
    fPropagators.push_back(prop);
    prop->fRunMgr = this;
    prop->SetConfig(fConfig);
    prop->fApplication      = fApplication;
    prop->fStdApplication   = fStdApplication;
    prop->fPhysicsInterface = fPhysicsInterface;
    prop->fPrimaryGenerator = fPrimaryGenerator;
    prop->fTruthMgr         = fTruthMgr;
  }

  // Temporary workaround to allow migration to detector construction
  if (fDetConstruction) {
    fDetConstruction->CreateMaterials();
    fDetConstruction->CreateGeometry();
    fNvolumes = fDetConstruction->SetupGeometry(fVolumes);
  } else {
    LoadGeometry(fConfig->fGeomFileName.c_str());
  }
  fConfig->fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
  Info(methodName, "Geometry created with maxdepth %d\n", fConfig->fMaxDepth);

  // Now we know the geometry depth: create the track data manager
  TrackDataMgr *dataMgr = TrackDataMgr::GetInstance(fConfig->fMaxDepth);

  // Initialize the process(es)
  if (!fPhysicsInterface) {
    geant::Fatal(methodName, "The physics process interface has to be initialized before this");
    return false;
  }
  // Initialize the physics
  fPhysicsInterface->Initialize();

  // Configure the locality manager and create the tracks
  LocalityManager *mgr = LocalityManager::Instance();
  if (!mgr->IsInitialized()) {
    mgr->SetNblocks(10);     // <- must be configurable
    mgr->SetBlockSize(1000); // <- must be configurable
    mgr->Init();
    if (fConfig->fUseVectorizedGeom) {
      printf("*** Using vectorized geometry, default basket size is %d\n", fConfig->fMaxPerBasket);
      if (fNthreads > 1) {
        printf("### \e[5mWARNING!    Basketized mode + MT not suported yet\e[m\n###\n");
      }
    } else {
      printf("*** Using scalar geometry\n");
    }
#if defined(GEANT_USE_NUMA) && !defined(VECCORE_CUDA_DEVICE_COMPILATION)
    if (fConfig->fUseNuma) {
      int nnodes = mgr->GetPolicy().GetNnumaNodes();
      mgr->SetPolicy(NumaPolicy::kCompact);
      // Loop propagatpors and assign NUMA nodes
      if (fNpropagators > 1) {
        for (auto i = 0; i < fNpropagators; ++i)
          fPropagators[i]->SetNuma(i % nnodes);
      }
    } else {
      mgr->SetPolicy(NumaPolicy::kSysDefault);
    }
#endif
  }

  fDoneEvents = BitSet::MakeInstance(fConfig->fNtotal);

  if (!fInitialisedRKIntegration) {
    PrepareRkIntegration();
  }

  if (fPrimaryGenerator) fPrimaryGenerator->InitPrimaryGenerator();

  // Initialize task data
  int nthreads = GetNthreadsTotal();
  fTDManager   = new TDManager(nthreads, fConfig->fMaxPerBasket);

  // Initialize application
  if (fConfig->fUseStdScoring) {
    fStdApplication = new StdApplication(this);
    fStdApplication->Initialize();
  }
  fApplication->Initialize();

  // Attach user data and physics data to task data
  for (int i = 0; i < nthreads; i++) {
    TaskData *td = fTDManager->GetTaskData(i);
    if (fPhysicsInterface) fPhysicsInterface->AttachUserData(td);
    if (fStdApplication) fStdApplication->AttachUserData(td);
    fApplication->AttachUserData(td);
  }

  // Initialize all propagators
  for (auto i = 0; i < fNpropagators; ++i)
    fPropagators[i]->Initialize();

  TaskData *td = fTDManager->GetTaskData(0);

  if (fConfig->fRunMode == GeantConfig::kExternalLoop) {
    for (auto i = 0; i < fNpropagators; ++i) {
      for (auto j = 0; j < fNthreads; ++j) {
        td = fTDManager->GetTaskData(i * fNthreads + j);
        td->AttachPropagator(fPropagators[i], 0);
      }
    }
  } else {
    td->AttachPropagator(fPropagators[0], 0);
  }

  // Initialize the event server
  fEventServer = new EventServer(fConfig->fNbuff, this);

  dataMgr->Print();
  fInitialized = true;
  return fInitialized;
}

//______________________________________________________________________________
TaskData *RunManager::BookTransportTask()
{
  // Book a transport task to be used with RunSimulationTask.
  return fTDManager->GetTaskData();
}

//______________________________________________________________________________
RunManager::~RunManager()
{
  for (auto i = 0; i < fNpropagators; ++i)
    delete fPropagators[i];
  fPropagators.clear();
  for (auto i = 0; i < fNvolumes; ++i) {
    Volume_t *vol          = (Volume_t *)fVolumes[i];
    VBconnector *connector = (VBconnector *)vol->GetBasketManagerPtr();
    vol->SetBasketManagerPtr(nullptr);
    delete connector;
  }
  BitSet::ReleaseInstance(fDoneEvents);
  delete fPhysicsInterface;

  // Cleanup user data attached to task data
  for (size_t i = 0; i < fTDManager->GetNtaskData(); ++i)
    fApplication->DeleteUserData(fTDManager->GetTaskData(i));
  delete fTDManager;
  delete fConfig;
  delete fStdApplication;
  delete fApplication;
}

//______________________________________________________________________________
bool RunManager::LoadGeometry(const char *filename)
{
// Load geometry from given file.
#ifdef USE_ROOT
  if (!gGeoManager) TGeoManager::Import(filename);
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    UserDetectorConstruction::LoadVecGeomGeometry(fBroker);
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
  } else {
    Error("Propagator::LoadGeometry", "Cannot load geometry from file %s", filename);
    return false;
  }
#else
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
  if (geom) {
    geom->LoadGeometryFromSharedLib(filename);
    vecgeom::GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
    fNvolumes = fVolumes.size();
  } else {
    Error("Propagator::LoadGeometry", "Cannot load geometry from file %s", filename);
    return false;
  }
#endif
  for (auto i = 0; i < fNvolumes; ++i) {
    Volume_t *vol          = (Volume_t *)fVolumes[i];
    VBconnector *connector = new VBconnector(i);
    vol->SetBasketManagerPtr(connector);
  }
  return true;
}

//______________________________________________________________________________
void RunManager::PrepareRkIntegration()
{
  using GUFieldPropagatorPool = ::GUFieldPropagatorPool;
  // using GUFieldPropagator = ::GUFieldPropagator;
  using std::cout;
  using std::endl;

  // Initialise the classes required for tracking in field
  // const unsigned int Nvar = 6; // Integration will occur over 3-position & 3-momentum coord.
  // const double hminimum = 1.0e-5; // * centimeter; =  0.0001 * millimeter;  // Minimum step = 0.1 microns
  // int statisticsVerbosity = 0;
  // cout << "Parameters for RK integration in magnetic field: " << endl;
  // cout << "   Driver parameters:  eps_tol= " << fConfig->fEpsilonRK << "  h_min= " << hminimum << endl;

  // auto integrDriver = new ScalarIntegrationDriver(hminimum, aStepper, Nvar, statisticsVerbosity);
  // auto fieldPropagator = new GUFieldPropagator(integrDriver, fConfig->fEpsilonRK);

  UserFieldConstruction *udc = fDetConstruction->GetFieldConstruction();
  if (!udc) {
    geant::Error("PrepareRkIntegration", "Cannot find expected User Field Construction object.");
    exit(1);
  } else {
    // GUVVectorField** fieldPtr;
    bool useRungeKutta = fConfig->fUseRungeKutta;
    cout << "PrepareRkIntegration: creating field & solver.  Use RK= " << useRungeKutta << endl;
    udc->CreateFieldAndSolver(useRungeKutta); // , fieldPtr );

    static GUFieldPropagatorPool *fpPool = GUFieldPropagatorPool::Instance();
    assert(fpPool); // Cannot be zero
    if (fpPool) {
      // fpPool->RegisterPrototype(fieldPropagator);
      //     ==> Now done in UserFieldConstructor::CreateSolverForField
      // Create clones for other threads
      fpPool->Initialize(fNthreads);
    } else {
      geant::Error("PrepareRkIntegration", "Cannot find GUFieldPropagatorPool Instance.");
    }
  }
}

//______________________________________________________________________________
void RunManager::SetUserFieldConstruction(UserFieldConstruction *udc)
{
  if (fDetConstruction)
    fDetConstruction->SetUserFieldConstruction(udc);
  else
    Error("RunManager::SetUserFieldConstruction",
          "To define a field, the user detector construction has to be defined");
  fInitialisedRKIntegration = false; //  Needs to be re-done !!
}

//______________________________________________________________________________
void RunManager::EventTransported(Event *event, TaskData *td)
{
  // Actions executed after an event is transported.
  // Adjust number of prioritized events
  if (event->IsPrioritized()) fPriorityEvents--;
  // event->Print();
  // Digitizer
  Info("EventTransported", " = task %d completed event %d with %d tracks", td->fTid, event->GetEvent(),
       event->GetNtracks());
  //  LocalityManager *lmgr = LocalityManager::Instance();
  //  Printf("   NQUEUED = %d  NBLOCKS = %d NRELEASED = %d",
  //         lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
  fApplication->FinishEvent(event);
  // Signal completion of one event to the event server
  fDoneEvents->SetBitNumber(event->GetEvent());
  assert(event->GetNtracks() > 0);
  // Notify event sets
  if (fConfig->fRunMode == GeantConfig::kExternalLoop) NotifyEventSets(event);

  int evtnb = event->GetEvent();

  fEventServer->CompletedEvent(event, td);

  // closing event in MCTruthManager
  if (fTruthMgr) fTruthMgr->CloseEvent(evtnb);
}

//______________________________________________________________________________
void RunManager::AddEventSet(EventSet *workload)
{
  // Add one event set to be processed in the system. Adds also the individual
  // events in the event server.

  LockEventSets();
  fEventSets.push_back(workload);
  UnlockEventSets();
  workload->AddSetToServer(fEventServer);
}

//______________________________________________________________________________
EventSet *RunManager::NotifyEventSets(Event *finished_event)
{
  // The method loops over registered event sets calling MarkDone method.
  LockEventSets();
  for (auto eventSet : fEventSets) {
    if (eventSet->MarkDone(finished_event->GetEvent())) {
      if (eventSet->IsDone()) {
        fEventSets.erase(std::remove(fEventSets.begin(), fEventSets.end(), eventSet), fEventSets.end());
        UnlockEventSets();
        return eventSet;
      }
      UnlockEventSets();
      return nullptr;
    }
  }
  UnlockEventSets();
  return nullptr;
}

//______________________________________________________________________________
bool RunManager::RunSimulationTask(EventSet *workload, TaskData *td)
{
  // Entry point for running simulation as asynchonous task. The user has to provide
  // an event set to be transported, and to pre-book the transport task. The method
  // will return only when the given workload is completed.
  //
  // The actual thread executing this task will co-operate
  // to the completion of other workloads pipelined by other similar tasks. In case
  // the worker doesn't manage to complete the workload and there is no more work
  // to be done, the thread will go to sleep and be waken at the completion of the
  // event set.
  // The method is intended to work in an external event loop scenario, steered by
  // an external task.

  // Register the workload in the manager and insert events in the server
  AddEventSet(workload);

  // LockEventSets();
  // std::cerr<<"Task "<< td->fTid <<" started event set:\n";
  // workload->Print();
  // UnlockEventSets();

  bool completed = WorkloadManager::TransportTracksTask(workload, td);

  // LockEventSets();
  // std::cerr<<"Task "<< td->fTid <<" finished event set:\n";
  // workload->Print();
  // UnlockEventSets();

  return completed;
}

//______________________________________________________________________________
void RunManager::RunSimulation()
{
  // Start simulation for all propagators
  Initialize();

  Printf("==========================================================================");
  Printf("= GeantV run started with %d propagator(s) using %d worker threads each ====", fNpropagators, fNthreads);
  if (fConfig->fUseRungeKutta)
    Printf("  Runge-Kutta integration ON with epsilon= %g", fConfig->fEpsilonRK);
  else
    Printf("  Runge-Kutta integration OFF");
  Printf("==========================================================================");
#ifdef USE_ROOT
  if (fConfig->fUseMonitoring) new TCanvas("cscheduler", "Scheduler monitor", 900, 600);
  if (fConfig->fUseAppMonitoring) new TCanvas("capp", "Application canvas", 700, 800);
#endif
  vecgeom::Stopwatch timer;
  timer.Start();
  for (auto i = 0; i < fNpropagators; ++i)
    fListThreads.emplace_back(Propagator::RunSimulation, fPropagators[i], fNthreads);

  for (auto &t : fListThreads) {
    t.join();
  }
  timer.Stop();
  double rtime      = timer.Elapsed();
  double ctime      = timer.CpuElapsed();
  long ntransported = fEventServer->GetNprimaries();
  long nsteps       = 0;
  long nsnext       = 0;
  long nphys        = 0;
  long nmag         = 0;
  long nsmall       = 0;
  long ncross       = 0;
  long npushed      = 0;
  long nkilled      = 0;

  for (auto i = 0; i < fNpropagators; ++i) {
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

  TaskData *td0 = fTDManager->GetTaskData(0);
  for (size_t stage = 0; stage < kNstages; ++stage) {
    SimulationStage *simstage = fPropagators[0]->GetStage(ESimulationStage(stage));
    if (!simstage->IsBasketized()) {
      Printf("Stage %20s: not basketized", simstage->GetName());
      continue;
    }
    // Merge stage counters
    for (size_t i = 1; i < fTDManager->GetNtaskData(); ++i) {
      TaskData *td = fTDManager->GetTaskData(i);
      *td0->fCounters[stage] += *td->fCounters[stage];
    }
    float nbasketized = td0->fCounters[stage]->fNvector;
    float ntotal      = nbasketized + td0->fCounters[stage]->fNscalar;
    Printf("Stage %20s: basketizable %ld/%d handlers > basketized %d %% (nscalar = %ld  nvector = %ld)",
           simstage->GetName(), simstage->GetNbasketized(), simstage->GetNhandlers(), int(100 * nbasketized / ntotal),
           size_t(ntotal - nbasketized), size_t(nbasketized));
  }

  Printf("=== Summary: %d propagators x %d threads: %ld primaries/%ld tracks,  total steps: %ld, snext calls: %ld, "
         "phys steps: %ld, mag. field steps: %ld, small steps: %ld, pushed: %ld, killed: %ld, bdr. crossings: %ld  "
         "RealTime=%gs CpuTime=%gs",
         fNpropagators, fNthreads, fEventServer->GetNprimaries(), ntransported, nsteps, nsnext, nphys, nmag, nsmall,
         npushed, nkilled, ncross, rtime, ctime);
  // LocalityManager *lmgr = LocalityManager::Instance();
  // Printf("NQUEUED = %d  NBLOCKS = %d NRELEASED = %d",
  //       lmgr->GetNqueued(), lmgr->GetNallocated(), lmgr->GetNreleased());
  FinishRun();
#ifdef USE_ROOT
  if (gApplication) delete gApplication;
#endif
}

//______________________________________________________________________________
bool RunManager::FinishRun()
{
  // Run termination actions.
  fApplication->FinishRun();
  if (fStdApplication) fStdApplication->FinishRun();
  // Actions to follow
  return true;
}

//______________________________________________________________________________
void RunManager::StopTransport()
{
  // Signal all propagators that transport has stopped
  for (auto i = 0; i < fNpropagators; ++i) {
    fPropagators[i]->StopTransport();
  }
}

} // GEANT_IMPL_NAMESPACE
} // Geant
