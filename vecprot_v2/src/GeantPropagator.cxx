// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The ProposeStep() method choses between a "scattering" process with no eloss
// and a "ionization" process and generates a random "physical" step. In this simple
// model all particlea undertake the same list of processes
// The ScatteringProcess() method emulates scattering and changes the particle
// direction randomly in a forward cone with opening angle proportional with 1/p
// The IonizationProcess() method simulates an energy deposition with an amount
// epsil*int(1+K*rnd) (epsil, 2*epsil, ..., K*epsil)
// In future we can add absorption and decay to check the stack model...
// The PropagateInField(step) method propagates in a uniform magnetic field with
// an amount equal to step
// The Transport() method computes the safety and snext and compares with the
// physical step. If (safety>pstep), PropagateInField(pstep) is called and the
// physics process is simulated. Otherwise, PropagateInField(safety) is called
// and safety subtracted from the physical step, then the procedure is repeated
// until C*snext/4 < 1E-6 (tangent of angle with sagita, C=1/R is the curvature)
//
#include "GeantPropagator.h"

#ifdef USE_ROOT
#include "TTimer.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <fenv.h>
#endif


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

#include "Geant/Error.h"
#include "GeantTrackVec.h"
#include "PhysicsInterface.h"
#include "PhysicsProcessOld.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "GeantVTaskMgr.h"
#include "StdApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"
#include "PrimaryGenerator.h"
#include "MCTruthMgr.h"

#ifdef USE_CALLGRIND_CONTROL
#include <valgrind/callgrind.h>
#endif

using namespace Geant;
using namespace vecgeom;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
    : fConfig(nullptr), fRunMgr(nullptr), fNthreads(1), fNbuff(100), fNtotal(1000), fNtransported(0),
      fNprimaries(0), fNsteps(0), fNsnext(0),
      fNphys(0), fNmag(0), fNsmall(0), fNcross(0),
      fTransportOngoing(false), fSingleTrack(false),
      fWMgr(nullptr), fApplication(nullptr), fStdApplication(nullptr),fTaskMgr(nullptr),
      fTimer(nullptr), fProcess(nullptr), fVectorPhysicsProcess(nullptr),
      fStoredTracks(nullptr), fPrimaryGenerator(nullptr), fTruthMgr(nullptr),
      fNtracks(0) {
  // Constructor 
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  delete fTimer;
  delete fWMgr;
}

//______________________________________________________________________________
int GeantPropagator::AddTrack(GeantTrack &track) {
  // Add a new track in the system. returns track number within the event.
  int slot = track.fEvslot;
  track.fParticle = fRunMgr->GetEvent(slot)->AddTrack();
  
  // call MCTruth manager if it has been instantiated
  if(fTruthMgr) fTruthMgr->AddTrack(track);

  //   fNtracks[slot]++;
  fNtransported++;
  return track.fParticle;
}

//______________________________________________________________________________
int GeantPropagator::DispatchTrack(GeantTrack &track, GeantTaskData *td) {
  // Dispatch a registered track produced by the generator.
  return fWMgr->GetScheduler()->AddTrack(track, td);
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(const GeantTrack_v &tracks, int itr) {
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);

  // stoping track in MCTruthManager
  if(fTruthMgr)
    {
      if(tracks.fStatusV[itr] == kKilled) fTruthMgr->EndTrack(tracks, itr);
    }
  
  if (fRunMgr->GetEvent(tracks.fEvslotV[itr])->StopTrack()) {
    std::atomic_int &priority_events = fRunMgr->GetPriorityEvents();
    priority_events++;
  }
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack(int tid) {
  // Returns a temporary track support for the physics processes, unique per
  // thread which can be used to add tracks produced by physics processes.
  if (tid < 0)
    tid = fWMgr->ThreadId();
  if (tid > fNthreads)
    Geant::Fatal("GeantPropagator::GetTempTrack", "Thread id %d is too large (max %d)", tid, fNthreads);
  GeantTrack &track = fRunMgr->GetTaskData(tid)->fTrack;
  track.Clear();
  return track;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantPropagator *GeantPropagator::NewInstance(int nthreads) {
  // Single instance of the propagator
  GeantPropagator* newInstance = new GeantPropagator();
  newInstance->fNthreads = nthreads;
  // Initialize workload manager
  #warning this must be solved before using CUDA
  newInstance->fWMgr = WorkloadManager::NewInstance(newInstance, nthreads);
  return newInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialize the propagator.
  // Add some empty baskets in the queue
  fWMgr->CreateBaskets(this);
}

//NOTE: We don't do anything here so it's not called from the WorkloadManager anymore
//______________________________________________________________________________
void GeantPropagator::ApplyMsc(int /*ntracks*/, GeantTrack_v & /*tracks*/, GeantTaskData * /*td*/) {
/*
  // Apply multiple scattering for charged particles.
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
#else
    mat = td->fVolume->GetMaterial();
#endif

#ifdef USE_REAL_PHYSICS
  tracks;
#else
  // actually nothing happens in this call in the TTabPhysProcess
  fProcess->ApplyMsc(mat, ntracks, tracks, td);
#endif
*/
}

//______________________________________________________________________________
void GeantPropagator::ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Generate all physics steps for the tracks in trackin.
  // Reset the current step length to 0
  for (int i = 0; i < ntracks; ++i) {
    tracks.fStepV[i] = 0.;
    tracks.fEdepV[i] = 0.;
  }
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = ((Medium_t *)td->fVolume->GetTrackingMediumPtr())->GetMaterial();
  ;
#else
    mat = td->fVolume->GetMaterial();
#endif

#ifdef USE_REAL_PHYSICS
  fPhysicsInterface->ComputeIntLen(mat, ntracks, tracks, 0, td);
#else
  fProcess->ComputeIntLen(mat, ntracks, tracks, 0, td);
#endif
}

//______________________________________________________________________________
void GeantPropagator::RunSimulation(GeantPropagator *prop, int nthreads, const char *geomfile, bool graphics, bool single)
{
// Static thread running the simulation for one propagator
  Info("RunSimulation", "Starting propagator %p with %d threads", prop, nthreads);
  prop->PropagatorGeom(geomfile, nthreads, graphics, single);
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, int nthreads, bool graphics, bool single) {
// Steering propagation method
  fConfig->fUseGraphics = graphics;
  fNthreads = nthreads;
  fSingleTrack = single;
  std::cout << " GeantPropagator::PropagatorGeom called with app= " << fApplication << std::endl;
  if (!fApplication) {
    printf("No user application attached - aborting");
    return;
  }

  fPrimaryGenerator->InitPrimaryGenerator();
  //   int itrack;

  if (fSingleTrack)
    Printf("==== Executing in single track loop mode using %d threads ====", fNthreads);
  else
    Printf("==== Executing in vectorized mode using %d threads ====", fNthreads);
  if (fConfig->fUsePhysics)
    Printf("  Physics ON");
  else
    Printf("  Physics OFF");
  if (fConfig->fUseRungeKutta)
    Printf("  Runge-Kutta integration ON with epsilon= %g", fConfig->fEpsilonRK);
  else
    Printf("  Runge-Kutta integration OFF");

  // Import the input events. This will start also populating the main queue
  if (!fEvents) {
    fEvents = new GeantEvent *[fNbuff];
    memset(fEvents, 0, fNbuff * sizeof(GeantEvent *));
  }

  // Loop baskets and transport particles until there is nothing to transport anymore
  fTransportOngoing = true;
  //fWMgr->SetMaxThreads(nthreads);
#ifdef USE_ROOT
  if (fConfig->fUseMonitoring) {
    TCanvas *cmon = new TCanvas("cscheduler", "Scheduler monitor", 900, 600);
    cmon->Update();
  }
  if (fConfig->fUseAppMonitoring) {
    TCanvas *capp = new TCanvas("capp", "Application canvas", 700, 800);
    capp->Update();
  }

  fTimer = new TStopwatch();
#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_START_INSTRUMENTATION;
#endif
#else
  fTimer = new vecgeom::Stopwatch();
#endif

  fTimer->Start();
  // Start system tasks
  if (!fWMgr->StartTasks(fTaskMgr)) {
    Fatal("PropagatorGeom", "Cannot start tasks.");
    return;
  }

  // Wake up the main scheduler once to avoid blocking the system
  //  condition_locker &sched_locker = fWMgr->GetSchLocker();
  //  sched_locker.StartOne();
  fWMgr->WaitWorkers();
  fWMgr->JoinThreads();
  fTimer->Stop();
#ifdef USE_CALLGRIND_CONTROL
  CALLGRIND_STOP_INSTRUMENTATION;
  CALLGRIND_DUMP_STATS;
#endif
#ifdef USE_ROOT
  double rtime = fTimer->RealTime();
  double ctime = fTimer->CpuTime();
#else
  double rtime = fTimer->Elapsed();
  double ctime = fTimer->Elapsed();  // TO FIX
#endif
  //   fTimer->Print();
  double speedup = ctime / rtime;
  double efficiency = speedup / nthreads;
//   fWMgr->Print();

#ifdef GEANTV_OUTPUT_RESULT_FILE
  const char *geomname = geomfile;
  if (strstr(geomfile, "http://root.cern.ch/files/"))
    geomname = geomfile + strlen("http://root.cern.ch/files/");
#endif
  //  int nsteps = fWMgr->GetScheduler()->GetNsteps();
  Printf("=== Transported: %ld primaries/%ld tracks,  total steps: %ld, snext calls: %ld, "
         "phys steps: %ld, mag. field steps: %ld, small steps: %ld bdr. crossings: %ld  RT=%gs, CP=%gs",
         fNprimaries.load(), fNtransported.load(), fNsteps.load(), fNsnext.load(), fNphys.load(), fNmag.load(),
         fNsmall.load(), fNcross.load(), rtime, ctime);
  Printf("   nthreads=%d speed-up=%f  efficiency=%f", nthreads, speedup, efficiency);
  //  Printf("Queue throughput: %g transactions/sec", double(fWMgr->FeederQueue()->n_ops()) / rtime);
  if (fTaskMgr) fTaskMgr->Finalize();
  fApplication->FinishRun();
  if (fStdApplication)
    fStdApplication->FinishRun();
#ifdef USE_VECGEOM_NAVIGATOR
  printf("=== Navigation done using VecGeom ====");
#else
  printf("=== Navigation done using TGeo    ====");
#endif
  //  Printf("Navstate pool usage statistics:");
  //   fWMgr->NavStates()->statistics();
}

//______________________________________________________________________________
void GeantPropagator::SetTaskBroker(TaskBroker *broker) {
  // Setter for task broker
  fWMgr->SetTaskBroker(broker);
}

//______________________________________________________________________________
TaskBroker *GeantPropagator::GetTaskBroker() {
  // Getter for task broker
  return fWMgr->GetTaskBroker();
}

//______________________________________________________________________________
void GeantPropagator::SetConfig(GeantConfig *config)
{
// Set run configuration.
  fConfig = config;
  fNbuff = config->fNbuff;
  fNtotal = config->fNtotal;
  if (fNtotal <= 0 || fNbuff <= 0) {
    Fatal("GeantPropagator::SetConfig", "Number of transported/buffered events should be positive");
    return;
  }
  if (fNbuff > fNtotal) {
    Info("GeantPropagator::SetCofig", "Number of buffered events changed to %d", fNtotal);
    fNbuff = fNtotal;
  }
  // Instantiate factory store
  GeantFactoryStore::Instance(fNbuff);
}
