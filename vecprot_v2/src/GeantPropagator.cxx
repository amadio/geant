// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The ProposeStep() method choses between a "scattering" process with no eloss
// and a "ionization" process and generates a random "physical" step. In this simple
// model all particlea undertake the same list of processes
// The ScatteringProcess() method emulates scattering and changes the particle
// direction randomly in a forward cone with opening angle proportional with 1/p
// The IonizationProcess() method simulates an energy deposition with an amount
// epsil*Int_t(1+K*rnd) (epsil, 2*epsil, ..., K*epsil)
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

#include "TSystem.h"
#include "TROOT.h"
#include "TTimer.h"
#include "TVirtualPad.h"
#include "TMath.h"
#include "TError.h"
#include "TGeoManager.h"
#include "TGeoHelix.h"
#include "TPolyMarker3D.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TGeoVolume.h"
#include "TGeoVoxelFinder.h"
#include "TGeoNode.h"
#include "TGeoMaterial.h"
#include <fenv.h>

#if USE_VECGEOM_NAVIGATOR == 1
#include "navigation/SimpleNavigator.h"
#include "management/RootGeoManager.h"
#include "volumes/LogicalVolume.h"
#include "volumes/PlacedVolume.h"
#endif
#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TGenPhaseSpace.h"
#include "GeantTrack.h"
#include "GeantOutput.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantThreadData.h"
#include "GeantVApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"
#include "PrimaryGenerator.h"

GeantPropagator *gPropagator = 0;

ClassImp(GeantPropagator)

    GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
    : TObject(), fNthreads(1), fNevents(100), fNtotal(1000), fNtransported(0), fNprimaries(0),
      fNsafeSteps(0), fNsnextSteps(0), fNphysSteps(0), fNprocesses(3), fNstart(0),
      fMaxTracks(0), // change
      fMaxThreads(100), fNminThreshold(10), fDebugTrk(-1), fMaxSteps(10000), fNperBasket(16),
      fMaxPerBasket(256), fMaxPerEvent(0), fMaxDepth(0), fLearnSteps(1000000), fMaxRes(10000.), fNaverage(0.), fVertex(),
      fEmin(1.E-4), // 100 KeV
      fEmax(10),    // 10 Gev
      fBmag(1.), fUsePhysics(kTRUE), fUseDebug(kFALSE), fUseGraphics(kFALSE),
      fTransportOngoing(kFALSE), fSingleTrack(kFALSE), fFillTree(kFALSE), fUseMonitoring(kFALSE),
      fUseAppMonitoring(kFALSE),
      fTracksLock(), fWMgr(0), fApplication(0), fOutput(0), fOutTree(0), fOutFile(0), fTimer(0),
      fProcess(0), fVectorPhysicsProcess(0), fStoredTracks(0), fPrimaryGenerator(0), fNtracks(0),
      fEvents(0), fThreadData(0) {
  // Constructor
  fVertex[0] = -8.;
  fVertex[1] = fVertex[2] = 0.;
  //   for (Int_t i=0; i<3; i++) fVertex[i] = gRandom->Gaus(0.,0.2);
  fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  Int_t i;
  delete fProcess;

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
  delete fOutput;
  delete fOutFile;
  delete fTimer;
  delete fWMgr;
  delete fApplication;
}

//______________________________________________________________________________
Int_t GeantPropagator::AddTrack(GeantTrack &track) {
  // Add a new track in the system. returns track number within the event.
  Int_t slot = track.fEvslot;
  track.fParticle = fEvents[slot]->AddTrack();
  //   fNtracks[slot]++;
  fNtransported++;
  return track.fParticle;
}

//______________________________________________________________________________
Int_t GeantPropagator::DispatchTrack(GeantTrack &track) {
  // Dispatch a registered track produced by the generator.
  return fWMgr->GetScheduler()->AddTrack(track);
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(const GeantTrack_v &tracks, Int_t itr) {
  // Mark track as stopped for tracking.
  //   Printf("Stopping track %d", track->particle);
  fEvents[tracks.fEvslotV[itr]]->StopTrack();
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack(Int_t tid) {
  // Returns a temporary track support for the physics processes, unique per
  // thread which can be used to add tracks produced by physics processes.
  if (tid < 0)
    tid = TGeoManager::ThreadId();
  if (tid > fNthreads)
    Fatal("GetTempTrack", "Thread id %d is too large (max %d)", tid, fNthreads);
  GeantTrack &track = fThreadData[tid]->fTrack;
  track.Clear();
  return track;
}

//______________________________________________________________________________
Int_t GeantPropagator::ImportTracks(Int_t nevents, Double_t average, Int_t startevent,
                                    Int_t startslot) {
  // Import tracks from "somewhere". Here we just generate nevents.
  static VolumePath_t *a = 0; // thread safe since initialized once used many times
#ifdef USE_VECGEOM_NAVIGATOR
  using vecgeom::SimpleNavigator;
  using vecgeom::Vector3D;
  using vecgeom::Precision;
  using vecgeom::GeoManager;
#endif

  Int_t tid = TGeoManager::ThreadId();
  if (tid > fNthreads)
    Fatal("ImportTracks", "Thread id %d is too large (max %d)", tid, fNthreads);
  GeantThreadData *td = fThreadData[tid];
  TGeoVolume *vol = 0;

  //   const Double_t etamin = -3, etamax = 3;
  Int_t ntracks = 0;
  Int_t ntotal = 0;
  Int_t ndispatched = 0;

  // the code below should be executed per track, as the primary vertex can change.
  if (!a) {
    a = VolumePath_t::MakeInstance(fMaxDepth);
#ifdef USE_VECGEOM_NAVIGATOR
    vecgeom::SimpleNavigator nav;
    nav.LocatePoint(GeoManager::Instance().GetWorld(),
                    Vector3D<Precision>(fVertex[0], fVertex[1], fVertex[2]), *a, true);
    vol = a->GetCurrentNode()->GetVolume();
    td->fVolume = vol;
#else
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    if (!nav)
      nav = gGeoManager->AddNavigator();
    TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
    // *td->fMatrix = nav->GetCurrentMatrix();
    vol = node->GetVolume();
    td->fVolume = vol;
    a->InitFromNavigator(nav);
#endif

  } else {
    TGeoNode const *node = a->GetCurrentNode();
    vol = node->GetVolume();
    td->fVolume = vol;
  }

  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
  Int_t threshold = nevents * average / (2 * fNthreads);
  threshold -= threshold % 4;
  if (threshold < 4)
    threshold = 4;
  if (threshold > fNperBasket)
    threshold = fNperBasket;
  basket_mgr->SetThreshold(threshold);

  static Bool_t init = kTRUE;
  if (init)
    init = kFALSE;
  Int_t event = startevent;
  for (Int_t slot = startslot; slot < startslot + nevents; slot++) {
    //     ntracks = td->fRndm->Poisson(average);
    ntracks = fPrimaryGenerator->NextEvent();

    ntotal += ntracks;
    fNprimaries += ntracks;
    if (!fEvents[slot])
      fEvents[slot] = new GeantEvent();
    fEvents[slot]->SetSlot(slot);
    fEvents[slot]->SetEvent(event);
    fEvents[slot]->Reset();

    for (Int_t i = 0; i < ntracks; i++) {
      GeantTrack &track = GetTempTrack(tid);
      track.SetPath(a);
      track.SetNextPath(a);
      track.SetEvent(event);
      track.SetEvslot(slot);
      fPrimaryGenerator->GetTrack(i, track);
      track.fFrombdr = kFALSE;
      track.fStatus = kAlive;
      track.fVindex = basket_mgr->GetNumber();
      AddTrack(track);
      ndispatched += DispatchTrack(track);
    }
    event++;
  }

  Printf("Imported %d tracks from events %d to %d. Dispatched %d baskets.", ntotal, startevent,
         startevent + nevents - 1, ndispatched);
  return ndispatched;
}

//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance(Int_t ntotal, Int_t nbuffered) {
  // Single instance of the propagator
  if (!fgInstance)
    fgInstance = new GeantPropagator();
  if (ntotal)
    fgInstance->fNtotal = ntotal;
  if (nbuffered) {
    fgInstance->fNevents = nbuffered;
    GeantFactoryStore::Instance(nbuffered);
  }
  return fgInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialization
  fMaxPerEvent = 5 * fNaverage;
  fMaxTracks = fMaxPerEvent * fNevents;

  // Initialize arrays here.
  gPropagator = GeantPropagator::Instance();

  if (!fProcess) {
    Fatal("Initialize", "The physics process has to be initialized before this");
    return;
  }
  // Initialize the process(es)
  fProcess->Initialize();
#if USE_VECPHYS == 1
  fVectorPhysicsProcess->Initialize();
#endif

  // Initialize workload manager
  fWMgr = WorkloadManager::Instance(fNthreads);
  if (fNthreads > fWMgr->GetNthreads()) {
    Error("Initialize", "Workload manager configured to support only %d thread but %d were "
                        "requested, using only %d threads.",
          fWMgr->GetNthreads(), fNthreads, fWMgr->GetNthreads());
    fNthreads = fWMgr->GetNthreads();
  }
  // Add some empty baskets in the queue
  fWMgr->CreateBaskets(); // geometry should be created by now

  if (!fNtracks) {
    fNtracks = new Int_t[fNevents];
    memset(fNtracks, 0, fNevents * sizeof(Int_t));
  }

  if (!fThreadData) {
    fThreadData = new GeantThreadData *[fNthreads + 1];
    for (Int_t i = 0; i < fNthreads + 1; i++)
      fThreadData[i] = new GeantThreadData(fMaxPerBasket, 3);
  }
  // Initialize application
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
    std::vector<vecgeom::LogicalVolume *> v1;
    vecgeom::GeoManager::Instance().getAllLogicalVolumes(v1);
    Printf("Have logical volumes %ld\n", v1.size());
    std::vector<vecgeom::VPlacedVolume *> v2;
    vecgeom::GeoManager::Instance().getAllPlacedVolumes(v2);
    Printf("Have placed volumes %ld\n", v2.size());
    vecgeom::RootGeoManager::Instance().world()->PrintContent();
  }
  return true;
}
#endif
//______________________________________________________________________________
Bool_t GeantPropagator::LoadGeometry(const char *filename) {
  // Load the detector geometry from file, unless already loaded.
  TGeoManager *geom = (gGeoManager) ? gGeoManager : TGeoManager::Import(filename);
  if (geom) {
#if USE_VECGEOM_NAVIGATOR == 1
    LoadVecGeomGeometry();
#endif
    fMaxDepth = TGeoManager::GetMaxLevels();
    return kTRUE;
  }
  ::Error("LoadGeometry", "Cannot load geometry from file %s", filename);
  return kFALSE;
}

//______________________________________________________________________________
void GeantPropagator::ApplyMsc(Int_t ntracks, GeantTrack_v &tracks, Int_t tid) {
  // Apply multiple scattering for charged particles.
  GeantThreadData *td = fThreadData[tid];
  TGeoMaterial *mat = 0;
  if (td->fVolume)
    mat = td->fVolume->GetMaterial();
  fProcess->ApplyMsc(mat, ntracks, tracks, tid);
}

//______________________________________________________________________________
void GeantPropagator::ProposeStep(Int_t ntracks, GeantTrack_v &tracks, Int_t tid) {
  // Generate all physics steps for the tracks in trackin.
  GeantThreadData *td = fThreadData[tid];
  // Reset the current step length to 0
  for (Int_t i = 0; i < ntracks; ++i) {
    tracks.fStepV[i] = 0.;
    tracks.fEdepV[i] = 0.;
  }
  TGeoMaterial *mat = 0;
  if (td->fVolume)
    mat = td->fVolume->GetMaterial();
  fProcess->ComputeIntLen(mat, ntracks, tracks, 0, tid);
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, Int_t nthreads, Bool_t graphics,
                                     Bool_t single) {
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
  // Initialize geometry and current volume
  if (!LoadGeometry(geomfile))
    return;
  Initialize();
  if (called) {
    Printf("Sorry, you can call this only once per session.");
    return;
  }
  called = kTRUE;

  fPrimaryGenerator->InitPrimaryGenerator();
  //   Int_t itrack;

  if (fSingleTrack)
    Printf("==== Executing in single track loop mode using %d threads ====", fNthreads);
  else
    Printf("==== Executing in vectorized mode using %d threads ====", fNthreads);
  if (fFillTree)
    Printf("  I/O enabled - disable if comparing single track loop with vectorized modes");
  else
    Printf("  I/O disabled");
  if (fUsePhysics)
    Printf("  Physics ON with %d processes", fNprocesses);
  else
    Printf("  Physics OFF");

  // Import the input events. This will start also populating the main queue
  if (!fEvents) {
    fEvents = new GeantEvent *[fNevents];
    memset(fEvents, 0, fNevents * sizeof(GeantEvent *));
  }

  ImportTracks(fNevents, fNaverage, 0, 0);

  // Initialize tree
  fOutput = new GeantOutput();
  fOutput->Init(fMaxTracks);
  if (fFillTree) {
    fOutFile = new TFile("output.root", "RECREATE");
    fOutTree = new TTree("TK", "Transport track data");
    fOutTree->Branch("gen", &fOutput);
  }

  // Loop baskets and transport particles until there is nothing to transport anymore
  fTransportOngoing = kTRUE;
  gGeoManager->SetMaxThreads(nthreads);
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
  condition_locker &sched_locker = fWMgr->GetSchLocker();
  sched_locker.StartOne();
  fWMgr->WaitWorkers();
  fTimer->Stop();
  Double_t rtime = fTimer->RealTime();
  Double_t ctime = fTimer->CpuTime();
  if (fFillTree)
    fOutTree->AutoSave();
  delete fOutFile;
  //   fTimer->Print();
  Double_t speedup = ctime / rtime;
  Double_t efficiency = speedup / nthreads;
  //   fWMgr->Print();
  fWMgr->JoinThreads();

#ifdef GEANTV_OUTPUT_RESULT_FILE
  const char *geomname = geomfile;
  if (strstr(geomfile, "http://root.cern.ch/files/"))
    geomname = geomfile + strlen("http://root.cern.ch/files/");
#endif
  Int_t nsteps = fWMgr->GetScheduler()->GetNsteps();
  Printf("=== Transported: %lld primaries/%lld tracks,  total steps: %d, safety steps: %lld,  snext steps: %lld, "
         "phys steps: %lld, RT=%gs, CP=%gs",
         fNprimaries.load(), fNtransported.load(), nsteps, fNsafeSteps.load(), fNsnextSteps.load(),
         fNphysSteps.load(), rtime, ctime);
  Printf("   nthreads=%d + 1 garbage collector speed-up=%f  efficiency=%f", nthreads, speedup,
         efficiency);
  Printf("Queue throughput: %g transactions/sec", double(fWMgr->FeederQueue()->n_ops()) / rtime);
#ifdef USE_VECGEOM_NAVIGATOR
  Printf("=== Navigation done using VecGeom ====");
#else
  Printf("=== Navigation done using TGeo    ====");
#endif
  Printf("Navstate pool usage statistics:");
//   fWMgr->NavStates()->statistics();
#ifdef GEANTV_OUTPUT_RESULT_FILE
  gSystem->mkdir("results");
  FILE *fp = fopen(Form("results/%s_%d.dat", geomname, single), "w");
  fprintf(fp, "%d %lld %lld %lld %g %g", single, fNsafeSteps, fNsnextSteps, fNphysSteps, rtime,
          ctime);
  fclose(fp);
#endif
  fOutFile = 0;
  fOutTree = 0;
}
