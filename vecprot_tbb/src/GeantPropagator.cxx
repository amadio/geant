// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The PhysicsSelect() method choses between a "scattering" process with no eloss
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

#include "GeantPropagator.h"

#include "TSystem.h"
#include "TROOT.h"
#include "TTimer.h"
#include "TVirtualPad.h"
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
#include "TGeoMatrix.h"
#include "TGeoBranchArray.h"
#include "TTree.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TBits.h"
#include "TGLSAViewer.h"
#include "TControlBar.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TGenPhaseSpace.h"
#include "GeantTrack.h"
#include "Event_Hit.h"
#include "GeantOutput.h"
#include "PhysicsProcess.h"
#include "GeantVolumeBasket.h"
#include "WorkloadManager.h"
#include "GeantThreadData.h"

#include "base/Global.h"
using vecgeom::kTwoPi;
using std::numeric_limits;

#include "CollDispTask.h"

#include "tbb/task_scheduler_init.h"
#include "tbb/tick_count.h"

void myObserver::on_scheduler_entry(bool is_worker) {
  Printf("-------------------------------------------------------------------------------------------Started");
}

void myObserver::on_scheduler_exit(bool is_worker) {
  Printf("-------------------------------------------------------------------------------------------Finished");
}

GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
    : fNthreads(1), fNevents(100), fNtotal(1000), fNtracks(0), fNadded(), fNtransported(), fNsafeSteps(0),
      fNsnextSteps(0), fNprocesses(3), fElossInd(0), fNstart(0), fMaxTracks(0), // change
      fMaxThreads(100), fDebugTrk(-1), fMaxSteps(10000), fNperBasket(10), fMaxPerBasket(100), fMaxPerEvent(0),
      fNaverage(0.), fEmin(0.1), // 100 MeV
      fEmax(10),                 // 10 Gev
      fBmag(1.), fUsePhysics(kTRUE), fUseDebug(kFALSE), fUseGraphics(kFALSE), fSingleTrack(kFALSE), fFillTree(kFALSE),
      fWMgr(0), fOutput(0), fKineTF1(0), fOutTree(0), fOutFile(0), fTimer(0), fProcesses(0), fTracks(0), fEvents(0),
      fNimportedEvents(), niter(),
      /*niter2(),
      niter3(),
      niter4(),*/
      pnTasksTotal(), ppTasksTotal(), dTasksTotal() {
  // Constructor
  gRandom->SetSeed();
  for (int i = 0; i < 3; i++)
    fVertex[i] = gRandom->Gaus(0., 10.);
  fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator() {
  // Destructor
  int i;
  if (fProcesses) {
    for (i = 0; i < fNprocesses; i++)
      delete fProcesses[i];
    delete[] fProcesses;
  }
  delete[] fTracks;

  if (fEvents) {
    for (i = 0; i < fNevents; i++)
      delete fEvents[i];
    delete[] fEvents;
  }

  delete fOutput;
  delete fOutFile;
  delete fTimer;
  delete fWMgr;
}

//______________________________________________________________________________
int GeantPropagator::AddTrack(GeantTrack *track) {
  // Add a new track in the system.

  fTracksLock.Lock();

  int slot = track->evslot;
  int itrack = fMaxPerEvent * slot + fNtracks[slot];
  track->particle = itrack;
  fEvents[slot]->AddTrack();
  fTracks[itrack] = track;
  fNtracks[slot]++;
  if (fNtracks[slot] == fMaxPerEvent) {
    Fatal("AddTrack", "No room to add track");
  }

  fTracksLock.UnLock();

  return itrack;
}

//______________________________________________________________________________
GeantTrack *GeantPropagator::AddTrack(int evslot) {
  // Reuse the next track booked for a given event slot.

  fTracksLock.Lock();

  int itrack = fMaxPerEvent * evslot + fNtracks[evslot];
  GeantTrack *track = fTracks[itrack];
  track->Reset();
  track->particle = itrack;
  fEvents[evslot]->AddTrack();
  fNtracks[evslot]++;
  //   printf("slot %d, track %d = %p\n", evslot, fNtracks[evslot], track);

  if (fNtracks[evslot] == fMaxPerEvent) {
    Fatal("AddTrack", "No room to add track");
  }

  fTracksLock.UnLock();

  return track;
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(GeantTrack *track) {
  // Mark track as stopped for tracking.
  //   Printf("Stop tr fr ev %d in slot %d. Left %d tr", fEvents[track->evslot]->event, track->evslot,
  //   fEvents[track->evslot]->TracksLeft());

  if (track->IsAlive()) {
    fEvents[track->evslot]->StopTrack();
    fNtransported++;
  } else {
    Printf("Stopping dead track.");
  }
  track->Kill();

  if (fEvents[track->evslot]->Transported()) {
    PerThread::reference TBBperThread = fTBBthreadData.local();
    TBBperThread.PushFinishedEvent(track->evslot);

    //      Printf ("Finished event %d in slot %d\tfNcurPriorEvents=%d", track->event, track->evslot, fNcurPriorEvents);
    Printf("Fin ev %d in sl %d", track->event, track->evslot);

    fEventsStatus[fEvents[track->evslot]->event]++; // In principle this can be moved into the PropTask
  }
}

//______________________________________________________________________________
int GeantPropagator::ImportTracks(int nevents, double average, int startevent, int startslot) {
  // Import tracks from "somewhere". Here we just generate nevents.
  PerThread::reference TBBperThread = fTBBthreadData.local();

  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
  TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
  *TBBperThread.fMatrix = nav->GetCurrentMatrix();
  TBBperThread.fVolume = node->GetVolume();
  GeantVolumeBasket *basket = (GeantVolumeBasket *)TBBperThread.fVolume->GetField();
  TGeoBranchArray a;
  a.InitFromNavigator(gGeoManager->GetCurrentNavigator());

  const double etamin = -3, etamax = 3;
  int ntracks = 0;
  int ntracksTotal = 0;

  // Species generated for the moment N, P, e, photon
  const int kMaxPart = 9;
  const int pdgGen[9] = {kPiPlus, kPiMinus, kProton, kProtonBar, kNeutron, kNeutronBar, kElectron, kPositron, kGamma};
  const double pdgRelProb[9] = {1., 1., 1., 1., 1., 1., 1., 1., 1.};
  const Species_t pdgSpec[9] = {kHadron, kHadron, kHadron, kHadron, kHadron, kHadron, kLepton, kLepton, kLepton};
  static double pdgProb[9] = {0.};
  int pdgCount[9] = {0};

  static Bool_t init = kTRUE;
  if (init) {
    pdgProb[0] = pdgRelProb[0];
    for (int i = 1; i < kMaxPart; ++i)
      pdgProb[i] = pdgProb[i - 1] + pdgRelProb[i];
    init = kFALSE;
  }
  int event = startevent;
  for (int slot = startslot; slot < startslot + nevents; slot++) {
    ntracks = TBBperThread.fRndm->Poisson(average);

    ntracksTotal += ntracks;

    Printf("Importing event %d into slot %d with %d tracks.", event, slot, ntracks);

    if (!fEvents[slot])
      fEvents[slot] = new GeantEvent();
    fEvents[slot]->slot = slot;
    fEvents[slot]->event = event;
    fEvents[slot]->Reset();
    fNtracks[slot] = 0;

    for (int i = 0; i < ntracks; i++) {
      GeantTrack *track = AddTrack(slot);
      *track->path = a;
      *track->nextpath = a;
      track->event = event;
      track->evslot = slot;
      double prob = TBBperThread.fRndm->Uniform(0., pdgProb[kMaxPart - 1]);
      track->pdg = 0;
      for (int j = 0; j < kMaxPart; ++j) {
        if (prob <= pdgProb[j]) {
          track->pdg = pdgGen[j];
          track->species = pdgSpec[j];
          //            Printf("Generating a %s",TDatabasePDG::Instance()->GetParticle(track->pdg)->GetName());
          pdgCount[j]++;
          break;
        }
      }
      if (!track->pdg)
        Fatal("ImportTracks", "No particle generated!");
      TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->pdg);
      track->charge = part->Charge() / 3.;
      track->mass = part->Mass();
      track->xpos = fVertex[0];
      track->ypos = fVertex[1];
      track->zpos = fVertex[2];
      track->e = fKineTF1->GetRandom() + track->mass;
      double p = sqrt((track->e - track->mass) * (track->e + track->mass));
      double eta = TBBperThread.fRndm->Uniform(etamin, etamax); // multiplicity is flat in rapidity
      double theta = 2 * atan(exp(-eta));
      // double theta = acos((1.-2.*gRandom->Rndm()));
      double phi = kTwoPi * TBBperThread.fRndm->Rndm();
      track->px = p * sin(theta) * cos(phi);
      track->py = p * sin(theta) * sin(phi);
      track->pz = p * cos(theta);
      track->frombdr = kFALSE;
      int itrack = track->particle;

      TBBperThread.fCollection->AddTrack(itrack, basket);

      fNstart++;
    }
    //      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
    event++;
    for (int i = 0; i < kMaxPart; i++) {
      //         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(),
      //         pdgCount[i]);
      pdgCount[i] = 0;
    }
  }

  fNadded += ntracksTotal;

  //   Printf("Injecting %d events...", nevents);
  InjectCollection(TBBperThread.fCollection);
  return ntracksTotal;
}

//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance() {
  // Single instance of the propagator
  if (!fgInstance)
    fgInstance = new GeantPropagator();
  return fgInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize() {
  // Initialization
  fMaxPerEvent = 5 * fNaverage;
  fMaxTracks = fMaxPerEvent * fNevents;

  // Initialize arrays here.
  if (!fKineTF1) {
    fKineTF1 = new TF1("fKineTF1", "gaus", fEmin, fEmax);
    fKineTF1->SetParameters(1, 3 * fEmin, 5);
  }

  if (!fProcesses) {
    fProcesses = new PhysicsProcess *[fNprocesses];
    fProcesses[0] = new ScatteringProcess("Scattering");
    fProcesses[1] = new ElossProcess("Eloss");
    fElossInd = 1;
    fProcesses[2] = new InteractionProcess("Interaction");
  }

  if (!fNtracks) {
    fNtracks = new int[fNevents];
    memset(fNtracks, 0, fNevents * sizeof(int));
  }
  if (!fTracks) {
    fTracks = new GeantTrack *[fMaxTracks];
    memset(fTracks, 0, fMaxTracks * sizeof(GeantTrack *));
    // Pre-book tracks (Thanks to olav.smorholm@cern.ch)
    fTracksStorage.reserve(fMaxTracks);
    for (int i = 0; i < fMaxTracks; i++) {
      fTracksStorage.push_back(GeantTrack(0));
      fTracks[i] = &(fTracksStorage.back());
    }
  }

  GeantThreadData td(fMaxPerBasket, 3);
  fTBBthreadData = PerThread(td);

  fDispThrPriority = 50;
  fEventsStatus = new tbb::atomic<int>[fNtotal];
  for (int evn = 0; evn < fNtotal; evn++)
    fEventsStatus[evn] = 0;
  fPriorityRange[0] = -1;
  fPriorityRange[1] = -1;

  niter = -1;

  fWMgr = WorkloadManager::Instance();
  // Add some empty baskets in the queue
  fWMgr->AddEmptyBaskets(1000);

  if (fUseGraphics) {
    // Int iter
    hnb = new TH1I("hnb", "number of baskets in the transport queue; iteration#", 500000, 0, 500000);
    htracks = new TH1I("htracks", "number of tracks/basket; iteration#", 500000, 0, 500000);
    numOfTracksTransportedInIter =
        new TH1I("numOfTracksTransportedInIter", "number of tracks transported; iteration#", 500000, 0, 500000);
    /*numOfnPropTasks = new TH1I("numOfnPropTasks", "number of normal PropTasks; iteration#", 500000, 0, 500000);
    numOfpPropTasks = new TH1I("numOfpPropTasks", "number of priority PropTasks; iteration#", 500000, 0, 500000);
    numOfDispTasks = new TH1I("numOfDispTasks", "number of CollDispTasks; iteration#", 500000, 0, 500000);*/

    // In time
    numOfTracksTransportedInTime =
        new TH1I("numOfTracksTransportedInTime", "Tracks transported; (1/5)second", 2000, 0, 2000);
    numOfPNtasks = new TH1I("numOfPNtasks", "Number of normal PropTasks running; (1/5)second", 2000, 0, 2000);
    numOfPPtasks = new TH1I("numOfPPtasks", "Number of priority PropTasks running; (1/5)second", 2000, 0, 2000);
    numOfDtasks = new TH1I("numOfDtasks", "Number of DispTasks running; (1/5)second", 2000, 0, 2000);
    sizeOfFQ = new TH1I("sizeOfFQ", "Size of feeder_queue; (1/5)second", 2000, 0, 2000);
    sizeOfPFQ = new TH1I("sizeOfPFQ", "Size of priority_feeder_queue; (1/5)second", 2000, 0, 2000);
    sizeOfCQ = new TH1I("sizeOfCQ", "Size of collector_queue; (1/5)second", 2000, 0, 2000);

    // Statistics
    numOfCollsPerTask =
        new TH1I("numOfCollsPerTask", "Number of CollDispTasks; Number of collections poped", 60, 0, 60);
    numOfTracksInBasket = new TH1I("numOfTracksInBasket", "Number of normal baskets; tracks per basket",
                                   fNperBasket * 2, 0, fNperBasket * 2);
    numOfTracksInPriorBasket = new TH1I("numOfTracksInPriorBasket", "Number of priority baskets; tracks per basket",
                                        fNperBasket * 2, 0, fNperBasket * 2);
    numOfTracksInColl =
        new TH1I("numOfTracksInColl", "Number of collection; tracks per collection", 2 * fNaverage, 0, 2 * fNaverage);

    hnb->SetFillColor(kRed);
    htracks->SetFillColor(kGreen);
    numOfTracksTransportedInIter->SetFillColor(kBlue);
    /*
          numOfnPropTasks->SetFillColor(kRed);
          numOfpPropTasks->SetFillColor(kGreen);
          numOfDispTasks->SetFillColor(kBlue);
    */
    numOfTracksTransportedInTime->SetFillColor(kAzure);
    numOfPNtasks->SetFillColor(kRed);
    numOfPPtasks->SetFillColor(kBlue);
    numOfDtasks->SetFillColor(kGreen);
    sizeOfFQ->SetFillColor(kGreen);
    sizeOfPFQ->SetFillColor(kMagenta);
    sizeOfCQ->SetFillColor(kYellow);

    numOfCollsPerTask->SetFillColor(kAzure);
    numOfTracksInBasket->SetFillColor(kCyan);
    numOfTracksInPriorBasket->SetFillColor(kOrange);
    numOfTracksInColl->SetFillColor(kPink);
  }
}

//______________________________________________________________________________
void GeantPropagator::InjectCollection(GeantTrackCollection *inColl) {
  // Inject collector handled by a single thread in the collector queue.

  GeantTrackCollection *newcoll;

  if (fWMgr->tbb_collector_empty_queue.try_pop(newcoll)) {
  } else
    newcoll = new GeantTrackCollection(100);

  *newcoll = *inColl;
  fWMgr->tbb_collector_queue.push(newcoll);
  inColl->Clear();
}

//______________________________________________________________________________
Bool_t GeantPropagator::LoadGeometry(const char *filename) {
  // Load the detector geometry from file.
  if (gGeoManager)
    return kTRUE;
  TGeoManager *geom = TGeoManager::Import(filename);
  if (geom) {
    // Create the basket array
    int nvols = gGeoManager->GetListOfVolumes()->GetEntries();
    fWMgr->CreateBaskets(2 * nvols);
    return kTRUE;
  }
  ::Error("LoadGeometry", "Cannot load geometry from file %s", filename);
  return kFALSE;
}

//______________________________________________________________________________
void GeantPropagator::PhysicsSelect(int ntracks, int *trackin) {
  // Generate all physics steps for the tracks in trackin.
  // Vectorized, except the unavoidable Sort()

  PerThread::reference TBBperThread = fTBBthreadData.local();

  static const double maxlen = numeric_limits<double>.max();
  double pstep;
  int ipart, iproc;
  GeantTrack *track;
  double *procStep;
  // Fill interaction lengths for all processes and all particles
  for (iproc = 0; iproc < fNprocesses; iproc++) {
    procStep = TBBperThread.GetProcStep(iproc);
    fProcesses[iproc]->ComputeIntLen(TBBperThread.fVolume, ntracks, trackin, procStep);
  }
  // Loop tracks and select process
  for (int i = 0; i < ntracks; i++) {
    ipart = trackin[i];
    track = fTracks[ipart];
    track->step = maxlen;
    track->process = -1;
    for (iproc = 0; iproc < fNprocesses; iproc++) {
      procStep = TBBperThread.GetProcStep(iproc);
      pstep = procStep[i];
      if (pstep < track->step) {
        track->step = pstep;
        track->process = iproc;
      }
    }
    if (fUseDebug && (fDebugTrk == ipart || fDebugTrk < 0)) {
      Printf("   PhysicsSelect: track #%d - process=%d pstep=%g", ipart, track->process, track->step);
      if (track->step > 1.E200) {
        Printf("xxx");
      }
    }
  }
}

//______________________________________________________________________________
void GeantPropagator::PrintParticles(int *trackin, int ntracks) {
  // Print the detailed particles list.
  Printf("================ THREAD : particles list");
  for (int i = 0; i < ntracks; i++) {
    fTracks[trackin[i]]->Print();
  }
}

//______________________________________________________________________________
void *GeantPropagator::GlobalObserver(void *arg) {
  Printf("======== Observer started ========");
  int recvSig;

  GeantPropagator *prop = GeantPropagator::Instance();
  int i = 0;

  while (1) {
    if (prop->observerSigsQueue.try_pop(recvSig)) {
      break;
    }
    gSystem->Sleep(200);
    prop->numOfPNtasks->Fill(i, prop->pnTasksRunning);
    prop->numOfPPtasks->Fill(i, prop->ppTasksRunning);
    prop->numOfDtasks->Fill(i, prop->dTasksRunning);
    prop->sizeOfFQ->Fill(i, prop->fWMgr->tbb_feeder_queue.size());
    prop->sizeOfPFQ->Fill(i, prop->fWMgr->tbb_feeder_priority_queue.size());
    prop->sizeOfCQ->Fill(i, prop->fWMgr->tbb_collector_queue.size());

    prop->numOfTracksTransportedInTime->Fill(i, prop->fNtransported);

    i++;
  }

  Printf("======== Observer finished ========");
  return 0;
}

//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, Bool_t graphics, Bool_t single) {
  // Propagate fNevents in the volume containing the vertex.
  // Simulate 2 physics processes up to exiting the current volume.
  static Bool_t called = kFALSE;
  fUseGraphics = graphics;
  fSingleTrack = single;
  Initialize();
  if (called) {
    Printf("Sorry, you can call this only once per session.");
    return;
  }
  called = kTRUE;

  // Initialize geometry and current volume
  if (!LoadGeometry(geomfile))
    return;
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
  // Create the basket array

  // Import the input events. This will start also populating the main queue
  if (!fEvents) {
    fEvents = new GeantEvent *[fNevents];
    memset(fEvents, 0, fNevents * sizeof(GeantEvent *));
  }

  int NinjectedTracks = ImportTracks(fNevents, fNaverage, 0, 0);
  fNimportedEvents += fNevents;

  // Initialize tree
  fOutput = new GeantOutput();
  fOutput->Init(fMaxTracks);
  if (fFillTree) {
    fOutFile = new TFile("output.root", "RECREATE");
    fOutTree = new TTree("TK", "Transport track data");
    fOutTree->Branch("gen", &fOutput);
  }

  // Loop baskets and transport particles until there is nothing to transport anymore
  gGeoManager->SetMaxThreads(fNthreads);
  tbb::task_scheduler_init init(fNthreads);

  fTimer = new TStopwatch();

  if (fUseGraphics) {
    TThread *observerThread = new TThread(GeantPropagator::GlobalObserver);
    observerThread->Run();
  }

  //   propObserver.observe();

  fTimer->Start();
  tbb::tick_count t0 = tick_count::now();

  CollDispTask &topTask = *new (task::allocate_root()) CollDispTask(1);
  task::spawn_root_and_wait(topTask);

  do {
    // Here we check if really all the events were transported
    int first_not_transported = fNtotal;
    for (int evn = 0; evn < fNtotal; evn++) {
      if (fEventsStatus[evn] == 0) {
        first_not_transported = evn;
        break;
      }
    }
    if (first_not_transported < fNtotal) {
      Printf("Not all events transported. First not transported: %d", first_not_transported);
      Printf("Feeder queue size: %d\t\tCollector queue size: %d", fWMgr->tbb_feeder_queue.size(),
             fWMgr->tbb_collector_queue.size());
      Printf("Collections waiting: %d\t\tTracks waiting: %d", fCollsWaiting, fTracksWaiting);
      Printf("feeder_queue size: %d", fWMgr->tbb_feeder_queue.size());
      Printf("feeder_priority_queue size: %d", fWMgr->tbb_feeder_priority_queue.size());

      fCollsWaiting = 0;
      fTracksWaiting = 0;

      fGarbageCollMode = kTRUE;
      CollDispTask &garbTask = *new (task::allocate_root()) CollDispTask(fWMgr->tbb_collector_queue.size());
      task::spawn_root_and_wait(garbTask);
    } else {
      Printf("ALL EVENTS FINISHED");
      break;
    }
  } while (1);

  fTimer->Stop();
  tbb::tick_count t1 = tick_count::now();

  if (fUseGraphics) {
    observerSigsQueue.push(1);
  }

  /*------------------------------------------------------------------------------------------------------*/
  // Just histograms filling
  ///*
  if (fUseGraphics) {
    TPad *pad1 = 0, *pad2 = 0, *pad3 = 0;
    TCanvas *c1 = 0;

    c1 = new TCanvas("c2", "c2", 1200, 1200);
    c1->Divide(1, 3);

    pad1 = (TPad *)c1->cd(1);
    hnb->SetStats(kFALSE);
    hnb->GetXaxis()->SetRangeUser(0, niter / 10);
    hnb->Draw();

    pad2 = (TPad *)c1->cd(2);
    htracks->SetStats(kFALSE);
    htracks->GetXaxis()->SetRangeUser(0, niter / 10);
    htracks->Draw();

    pad3 = (TPad *)c1->cd(3);
    numOfTracksTransportedInIter->SetStats(kFALSE);
    numOfTracksTransportedInIter->GetXaxis()->SetRangeUser(0, niter / 10);
    numOfTracksTransportedInIter->Draw();

    TFile *diagFile = new TFile("diagnostics.root", "UPDATE");

    hnb->Write();
    htracks->Write();
    numOfTracksTransportedInIter->Write();
    /*
          numOfnPropTasks->Write();
          numOfpPropTasks->Write();
          numOfDispTasks->Write();
    */
    numOfPNtasks->Write();
    numOfPPtasks->Write();
    numOfDtasks->Write();
    sizeOfFQ->Write();
    sizeOfPFQ->Write();
    sizeOfCQ->Write();
    numOfCollsPerTask->Write();
    numOfTracksInBasket->Write();
    numOfTracksInPriorBasket->Write();
    numOfTracksInColl->Write();
    numOfTracksTransportedInTime->Write();

    diagFile->Close();
  }
  //*/
  /*------------------------------------------------------------------------------------------------------*/

  double rtime = fTimer->RealTime();
  double ctime = fTimer->CpuTime();
  double tbbtime = (t1 - t0).seconds();

  if (fFillTree)
    fOutTree->AutoSave();
  delete fOutFile;

  FILE *fp = fopen("rez.txt", "a");
  fprintf(fp, "geometry \"%s\"\n", geomfile);
  fprintf(fp, "threads %d\nevtotal %d\nevbuf %d\ntr/ev %f\ntr/bask %d\n", fNthreads, fNtotal, fNevents, fNaverage,
          fNperBasket);
  fprintf(fp, "added %lld\ntransported %lld\nsafesteps %lld\nsnextsteps %lld\nrealtime %g\ncputime %g\ntbbtime %g\n",
          fNadded, fNtransported, fNsafeSteps, fNsnextSteps, rtime, ctime, tbbtime);
  fprintf(fp, "Total %d normal PropTasks\nTotal %d priority PropTasks\nTotal %d DispTasks\n", pnTasksTotal,
          ppTasksTotal, dTasksTotal);
  fprintf(fp,
          "------------------------------------------------------------------------------------------------------\n");
  fclose(fp);

  FILE *fp2 = fopen("diag.txt", "a");
  fprintf(fp2, "%f, %g, %g, %f, %f\n", (double)fNthreads, ctime, rtime, (double)fNsafeSteps,
          (double)fNsnextSteps);
  fclose(fp2);

  fOutFile = 0;
  fOutTree = 0;
}

//______________________________________________________________________________
void GeantPropagator::SelectTracksForProcess(int iproc, int ntotransport, int *particles, int &ntodo,
                                             int *parttodo) {
  // Add to output array all particles that will do the process iproc.
  for (int itr = 0; itr < ntotransport; itr++)
    if (fTracks[particles[itr]]->process == iproc)
      parttodo[ntodo++] = particles[itr];
}
