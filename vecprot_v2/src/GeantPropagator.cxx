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

#include "TTimer.h"
#include "TError.h"
#include "TStopwatch.h"
#include "TCanvas.h"
#include <fenv.h>

#if USE_VECGEOM_NAVIGATOR == 1
#include "navigation/SimpleNavigator.h"
#include "management/RootGeoManager.h"
#include "volumes/PlacedVolume.h"
#else
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#endif

#include "GeantTrack.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantTaskData.h"
#include "GeantVApplication.h"
#include "StdApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"
#include "PrimaryGenerator.h"

using namespace Geant;

GeantPropagator *gPropagator = 0;

ClassImp(GeantPropagator)

    GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
    : TObject(), fNthreads(1), fNevents(100), fNtotal(1000), fNtransported(0), fNprimaries(0), fNsafeSteps(0),
      fNsnextSteps(0), fNphysSteps(0), fFeederLock(ATOMIC_FLAG_INIT), fPriorityEvents(0), fDoneEvents(0),
      fNprocesses(3), fNstart(0), fMaxTracks(0), fMaxThreads(100), fNminThreshold(10), fDebugEvt(-1), fDebugTrk(-1),
      fDebugStp(-1), fDebugRep(-1), fMaxSteps(10000), fNperBasket(16), fMaxPerBasket(256), fMaxPerEvent(0),
      fMaxDepth(0), fLearnSteps(0), fLastEvent(0), fPriorityThr(0), fMaxRes(0), fNaverage(0), fVertex(),
      fEmin(1.E-4), // 100 KeV
      fEmax(10),    // 10 Gev
      fBmag(1.), fUsePhysics(kTRUE), fUseDebug(kFALSE), fUseGraphics(kFALSE), fUseStdScoring(kFALSE),
      fTransportOngoing(kFALSE), fSingleTrack(kFALSE), fUseMonitoring(kFALSE), fUseAppMonitoring(kFALSE), fTracksLock(),
      fWMgr(0), fApplication(0), fStdApplication(0), fTimer(0), fProcess(0), fVectorPhysicsProcess(0), fStoredTracks(0),
      fPrimaryGenerator(0), fNtracks(0), fEvents(0), fThreadData(0) {
  // Constructor
  fVertex[0] = fVertex[1] = fVertex[2] = 0.;
  fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  int i;
  delete fProcess;
  BitSet::ReleaseInstance(fDoneEvents);
#if USE_VECPHYS == 1
  delete fVectorPhysicsProcess;
#endif

  if (fEvents) {
    for (i = 0; i < fNevents; i++)
      delete fEvents[i];
    delete[] fEvents;
  }

  if (fThreadData) {
    for (i = 0; i < fNthreads; i++)
      delete fThreadData[i];
    delete[] fThreadData;
  }
  delete fTimer;
  delete fWMgr;
  delete fApplication;
}

//______________________________________________________________________________
int GeantPropagator::AddTrack(GeantTrack &track) {
  // Add a new track in the system. returns track number within the event.
  int slot = track.fEvslot;
  track.fParticle = fEvents[slot]->AddTrack();
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
  if (fEvents[tracks.fEvslotV[itr]]->StopTrack())
    fPriorityEvents++;
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack(int tid) {
  // Returns a temporary track support for the physics processes, unique per
  // thread which can be used to add tracks produced by physics processes.
  if (tid < 0)
    tid = WorkloadManager::Instance()->ThreadId();
  if (tid > fNthreads)
    Fatal("GetTempTrack", "Thread id %d is too large (max %d)", tid, fNthreads);
  GeantTrack &track = fThreadData[tid]->fTrack;
  track.Clear();
  return track;
}

//______________________________________________________________________________
int GeantPropagator::Feeder(GeantTaskData *td) {
  // Feeder called by any thread to inject the next event(s)
  // Only one thread at a time
  if (fFeederLock.test_and_set(std::memory_order_acquire))
    return -1;
  int nbaskets = 0;
  if (!fLastEvent) {
    nbaskets = ImportTracks(fNevents, 0, 0, td);
    fLastEvent = fNevents;
    fFeederLock.clear(std::memory_order_release);
    return nbaskets;
  }
  // Check and mark finished events
  for (int islot = 0; islot < fNevents; islot++) {
    GeantEvent *evt = fEvents[islot];
    if (fDoneEvents->TestBitNumber(evt->GetEvent()))
      continue;
    if (evt->Transported()) {
      fPriorityEvents--;
      evt->Print();
      // Digitizer (todo)
      int ntracks = fNtracks[islot];
      Printf("= digitizing event %d with %d tracks pri=%d", evt->GetEvent(), ntracks, fPriorityEvents.load());
      //            propagator->fApplication->Digitize(evt->GetEvent());
      fDoneEvents->SetBitNumber(evt->GetEvent());
      if (fLastEvent < fNtotal) {
        Printf("=> Importing event %d", fLastEvent);
        nbaskets += ImportTracks(1, fLastEvent, islot, td);
        fLastEvent++;
      }
    }
  }

  fFeederLock.clear(std::memory_order_release);
  return nbaskets;
}

//______________________________________________________________________________
int GeantPropagator::ImportTracks(int nevents, int startevent, int startslot, GeantTaskData *thread_data) {
  // Import tracks from "somewhere". Here we just generate nevents.
  static VolumePath_t *a = 0; // thread safe since initialized once used many times
#ifdef USE_VECGEOM_NAVIGATOR
  using vecgeom::SimpleNavigator;
  using vecgeom::Vector3D;
  using vecgeom::Precision;
  using vecgeom::GeoManager;
#endif

  Volume_t *vol = 0;
  int ntracks = 0;
  int ntotal = 0;
  int ndispatched = 0;
  GeantTaskData *td = thread_data;
  if (td == 0) {
    int tid = WorkloadManager::Instance()->ThreadId();
    td = fThreadData[tid];
    td->fTid = tid;
  }
  // the code below should be executed per track, as the primary vertex can change.
  if (!a) {
    a = VolumePath_t::MakeInstance(fMaxDepth);
#ifdef USE_VECGEOM_NAVIGATOR
    vecgeom::SimpleNavigator nav;
    nav.LocatePoint(GeoManager::Instance().GetWorld(), Vector3D<Precision>(fVertex[0], fVertex[1], fVertex[2]), *a,
                    true);
    vol = const_cast<Volume_t *>(a->Top()->GetLogicalVolume());
    td->fVolume = vol;
#else
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    if (!nav)
      nav = gGeoManager->AddNavigator();
    TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
    vol = node->GetVolume();
    td->fVolume = vol;
    a->InitFromNavigator(nav);
#endif

  } else {
#ifdef USE_VECGEOM_NAVIGATOR
    vol = const_cast<Volume_t *>(a->Top()->GetLogicalVolume());
#else
    TGeoNode const *node = a->GetCurrentNode();
    vol = node->GetVolume();
#endif
    td->fVolume = vol;
  }

#ifdef USE_VECGEOM_NAVIGATOR
  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->getBasketManagerPtr());
#else
  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
#endif
  basket_mgr->SetThreshold(fNperBasket);

  static Bool_t init = kTRUE;
  if (init)
    init = kFALSE;
  int event = startevent;
  for (int slot = startslot; slot < startslot + nevents; slot++) {
    ntracks = fPrimaryGenerator->NextEvent();
    ntotal += ntracks;
    fNprimaries += ntracks;
    if (!fEvents[slot])
      fEvents[slot] = new GeantEvent();
    fEvents[slot]->SetSlot(slot);
    fEvents[slot]->SetEvent(event);
    fEvents[slot]->Reset();
    // Set priority threshold to non-default value
    if (fPriorityThr > 0)
      fEvents[slot]->SetPriorityThr(fPriorityThr);

    for (int i = 0; i < ntracks; i++) {
      GeantTrack &track = td->GetTrack();
      track.SetPath(a);
      track.SetNextPath(a);
      track.SetEvent(event);
      track.SetEvslot(slot);
      fPrimaryGenerator->GetTrack(i, track);
      if (!track.IsNormalized())
        track.Print();
      track.fFrombdr = kFALSE;
      track.fStatus = kAlive;
      track.fVindex = basket_mgr->GetNumber();
      AddTrack(track);
      ndispatched += DispatchTrack(track, thread_data);
    }
    event++;
  }

  Printf("Imported %d tracks from events %d to %d. Dispatched %d baskets.", ntotal, startevent,
         startevent + nevents - 1, ndispatched);
  return ndispatched;
}

//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance(int ntotal, int nbuffered, int nthreads) {
  // Single instance of the propagator
  if (fgInstance)
    return fgInstance;
  if (ntotal <= 0 || nbuffered <= 0) {
    Printf("GeantPropagator::Instance: Number of transported/buffered events should be positive");
    return 0;
  }
  fgInstance = new GeantPropagator();
  fgInstance->fNtotal = ntotal;
  fgInstance->fNevents = nbuffered;
  fgInstance->fNthreads = nthreads;
  if (nbuffered > ntotal) {
    Printf("GeantPropagator::Instance: Number of buffered events changed to %d", ntotal);
    fgInstance->fNevents = ntotal;
  }
  // Initialize workload manager
  fgInstance->fWMgr = WorkloadManager::Instance(nthreads);
  // Instantiate factory store
  GeantFactoryStore::Instance(nbuffered);
  return fgInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialization
  fMaxPerEvent = 5 * fNaverage;
  fMaxTracks = fMaxPerEvent * fNevents;

  // Initialize arrays here.
  gPropagator = GeantPropagator::Instance();
  fDoneEvents = BitSet::MakeInstance(fNtotal);
  if (!fProcess) {
    Fatal("Initialize", "The physics process has to be initialized before this");
    return;
  }
  // Initialize the process(es)
  fProcess->Initialize();
#if USE_VECPHYS == 1
  fVectorPhysicsProcess->Initialize();
#endif

  if (!fNtracks) {
    fNtracks = new int[fNevents];
    memset(fNtracks, 0, fNevents * sizeof(int));
  }
}

//______________________________________________________________________________
void GeantPropagator::InitializeAfterGeom() {
  // Initialization, part two.

  // Add some empty baskets in the queue
  fWMgr->CreateBaskets(); // geometry must be created by now

  if (!fThreadData) {
    fThreadData = new GeantTaskData *[fNthreads];
    for (int i = 0; i < fNthreads; i++) {
      fThreadData[i] = new GeantTaskData();
      fThreadData[i]->fTid = i;
      //      for (int j=0; j<1000; j++) {
      //        GeantBasket *b = new GeantBasket(fNperBasket, fMaxDepth);
      //        fThreadData[i]->RecycleBasket(b);
      //      }
    }
  }
  // Initialize application
  if (fUseStdScoring) {
    fStdApplication = new StdApplication();
    fStdApplication->Initialize();
  }
  fApplication->Initialize();
}

#if USE_VECGEOM_NAVIGATOR == 1
/**
 * function to setup the VecGeom geometry from a TGeo geometry ( if gGeoManager ) exists
 */
//______________________________________________________________________________
Bool_t GeantPropagator::LoadVecGeomGeometry() {
  if (vecgeom::GeoManager::Instance().GetWorld() == NULL) {
    Printf("Now loading VecGeom geometry\n");
    vecgeom::RootGeoManager::Instance().LoadRootGeometry();
    Printf("Loading VecGeom geometry done\n");
    Printf("Have depth %d\n", vecgeom::GeoManager::Instance().getMaxDepth());
    // Create the tab. phys process.
    Printf("Have logical volumes %d\n", vecgeom::GeoManager::Instance().GetLogicalVolumesCount());
    Printf("Have placed volumes %d\n", vecgeom::GeoManager::Instance().GetPlacedVolumesCount());
#ifdef USE_VECGEOM_NAVIGATOR
    std::cout << vecgeom::GeoManager::Instance().GetWorld() << endl;
#else
    vecgeom::RootGeoManager::Instance().world()->PrintContent();
#endif
  }
  if (fWMgr && fWMgr->GetTaskBroker()) {
    Printf("Now upload VecGeom geometry to Coprocessor(s)\n");
    return fWMgr->LoadGeometry();
  } else
    return true;
}
#endif

#ifndef USE_VECGEOM_GEOMETRY
//______________________________________________________________________________
Bool_t GeantPropagator::LoadGeometry(const char *filename) {
// Load the detector geometry from file, unless already loaded.
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::GeoManager *geom = &vecgeom::GeoManager::Instance();
#else
  TGeoManager *geom = (gGeoManager) ? gGeoManager : TGeoManager::Import(filename);
#endif
  if (geom) {
#if USE_VECGEOM_NAVIGATOR == 1
    LoadVecGeomGeometry();
    fMaxDepth = vecgeom::GeoManager::Instance().getMaxDepth();
#else
    fMaxDepth = TGeoManager::GetMaxLevels();
#endif
    return kTRUE;
  }
  ::Error("LoadGeometry", "Cannot load geometry from file %s", filename);
  return kFALSE;
}
#endif

//______________________________________________________________________________
void GeantPropagator::ApplyMsc(int ntracks, GeantTrack_v &tracks, GeantTaskData *td) {
  // Apply multiple scattering for charged particles.
  Material_t *mat = 0;
  if (td->fVolume)
#ifdef USE_VECGEOM_NAVIGATOR
    mat = ((Medium_t *)td->fVolume->getTrackingMediumPtr())->GetMaterial();
#else
    mat = td->fVolume->GetMaterial();
#endif
  fProcess->ApplyMsc(mat, ntracks, tracks, td);
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
    mat = ((Medium_t *)td->fVolume->getTrackingMediumPtr())->GetMaterial();
  ;
#else
    mat = td->fVolume->GetMaterial();
#endif
  fProcess->ComputeIntLen(mat, ntracks, tracks, 0, td);
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, int nthreads, Bool_t graphics, Bool_t single) {
  // Propagate fNevents in the volume containing the vertex.
  // Simulate 2 physics processes up to exiting the current volume.
  static Bool_t called = kFALSE;
  fUseGraphics = graphics;
  fNthreads = nthreads;
  fSingleTrack = single;
  if (!fApplication) {
    Printf("No user application attached - aborting");
    return;
  }
  Initialize();
  // Initialize geometry and current volume
  if (!LoadGeometry(geomfile))
    return;
  InitializeAfterGeom();
  // Initialize application
  fApplication->Initialize();
  if (called) {
    Printf("Sorry, you can call this only once per session.");
    return;
  }
  called = kTRUE;

  fPrimaryGenerator->InitPrimaryGenerator();
  //   int itrack;

  if (fSingleTrack)
    Printf("==== Executing in single track loop mode using %d threads ====", fNthreads);
  else
    Printf("==== Executing in vectorized mode using %d threads ====", fNthreads);
  if (fUsePhysics)
    Printf("  Physics ON with %d processes", fNprocesses);
  else
    Printf("  Physics OFF");

  // Import the input events. This will start also populating the main queue
  if (!fEvents) {
    fEvents = new GeantEvent *[fNevents];
    memset(fEvents, 0, fNevents * sizeof(GeantEvent *));
  }

  //  Feeder(fThreadData[0]);

  // Loop baskets and transport particles until there is nothing to transport anymore
  fTransportOngoing = kTRUE;
  WorkloadManager::Instance()->SetMaxThreads(nthreads);
  if (fUseMonitoring) {
    TCanvas *cmon = new TCanvas("cscheduler", "Scheduler monitor", 900, 600);
    cmon->Update();
  }
  if (fUseAppMonitoring) {
    TCanvas *capp = new TCanvas("capp", "Application canvas", 700, 800);
    capp->Update();
  }
  fTimer = new TStopwatch();
  fWMgr->StartThreads();
  fTimer->Start();
  // Wake up the main scheduler once to avoid blocking the system
  //  condition_locker &sched_locker = fWMgr->GetSchLocker();
  //  sched_locker.StartOne();
  fWMgr->WaitWorkers();
  fTimer->Stop();
  double rtime = fTimer->RealTime();
  double ctime = fTimer->CpuTime();
  //   fTimer->Print();
  double speedup = ctime / rtime;
  double efficiency = speedup / nthreads;
  //   fWMgr->Print();
  fWMgr->JoinThreads();

#ifdef GEANTV_OUTPUT_RESULT_FILE
  const char *geomname = geomfile;
  if (strstr(geomfile, "http://root.cern.ch/files/"))
    geomname = geomfile + strlen("http://root.cern.ch/files/");
#endif
  int nsteps = fWMgr->GetScheduler()->GetNsteps();
  Printf("=== Transported: %lld primaries/%lld tracks,  total steps: %d, safety steps: %lld,  snext steps: %lld, "
         "phys steps: %lld, RT=%gs, CP=%gs",
         fNprimaries.load(), fNtransported.load(), nsteps, fNsafeSteps.load(), fNsnextSteps.load(), fNphysSteps.load(),
         rtime, ctime);
  Printf("   nthreads=%d speed-up=%f  efficiency=%f", nthreads, speedup, efficiency);
  //  Printf("Queue throughput: %g transactions/sec", double(fWMgr->FeederQueue()->n_ops()) / rtime);
  fApplication->FinishRun();
  if (fStdApplication)
    fStdApplication->FinishRun();
#ifdef USE_VECGEOM_NAVIGATOR
  Printf("=== Navigation done using VecGeom ====");
#else
  Printf("=== Navigation done using TGeo    ====");
#endif
  //  Printf("Navstate pool usage statistics:");
  //   fWMgr->NavStates()->statistics();
}

//______________________________________________________________________________
int GeantPropagator::GetMonFeatures() const {
  // Get the number of monitored features
  return fWMgr->GetMonFeatures();
}

//______________________________________________________________________________
bool GeantPropagator::IsMonitored(EGeantMonitoringType feature) const {
  // Check if a given feature is monitored
  return fWMgr->IsMonitored(feature);
}

//______________________________________________________________________________
void GeantPropagator::SetMonitored(EGeantMonitoringType feature, bool flag) {
  // Enable monitoring a feature
  fWMgr->SetMonitored(feature, flag);
}

//______________________________________________________________________________
void GeantPropagator::SetNminThreshold(int thr) {
  // Setter for the global transport threshold
  fWMgr->SetNminThreshold(thr);
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
