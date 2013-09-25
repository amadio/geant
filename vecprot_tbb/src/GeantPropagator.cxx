// A simple propagator taking as input a set of particles located in a given
// volume AND the global matrix of the volume.
// The PhysicsSelect() method choses between a "scattering" process with no eloss
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

#include "CollDispTask.h"

#include "tbb/task_scheduler_init.h"

GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator() :
                 fNthreads(1),
                 fNevents(100),
                 fNtotal(1000),
                 fNtracks(0),
                 fNadded(),
                 fNtransported(),
                 fNsafeSteps(0),
                 fNsnextSteps(0),
                 fNprocesses(3),
                 fElossInd(0),
                 fNstart(0),
                 fMaxTracks(0),  //change
                 fMaxThreads(100),
                 fDebugTrk(-1),
                 fMaxSteps(10000),
                 fNperBasket(10),
                 fMaxPerBasket(100),
                 fMaxPerEvent(0),
                 fNaverage(0.),
                 fEmin(0.1), // 100 MeV
                 fEmax(10),  // 10 Gev
                 fBmag(1.),
                 fUsePhysics(kTRUE),
                 fUseDebug(kFALSE),
                 fUseGraphics(kFALSE),
                 fSingleTrack(kFALSE),
                 fFillTree(kFALSE),
                 fWMgr(0),
                 fOutput(0),
                 fKineTF1(0),
                 fOutTree(0),
                 fOutFile(0),
                 fTimer(0),
                 fProcesses(0),
                 fTracks(0),
                 fEvents(0),
                 fPrioritize(),
                 fNimportedEvents(),
                 niter()
{
// Constructor
   for (Int_t i=0; i<3; i++) fVertex[i] = gRandom->Gaus(0.,10.);
   fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator()
{
// Destructor
   Int_t i;
   if (fProcesses) {
     for (i=0; i<fNprocesses; i++) delete fProcesses[i];
     delete [] fProcesses;
   }
   delete [] fTracks;

   if (fEvents) {
      for (i=0; i<fNevents; i++) delete fEvents[i];
      delete [] fEvents;
   }

   delete fOutput;
   delete fOutFile;
   delete fTimer;
   delete fWMgr;
}

//______________________________________________________________________________
Int_t GeantPropagator::AddTrack(GeantTrack *track)
{
// Add a new track in the system.

   fTracksLock.Lock();

   Int_t slot = track->evslot;
   Int_t itrack = fMaxPerEvent*slot+fNtracks[slot];
   track->particle = itrack;
   fEvents[slot]->AddTrack();
   fTracks[itrack] = track;
   fNtracks[slot]++;
   if (fNtracks[slot]==fMaxPerEvent) {
      Fatal("AddTrack", "No room to add track");
   }

   fTracksLock.UnLock();

   return itrack;
}

//______________________________________________________________________________
GeantTrack *GeantPropagator::AddTrack(Int_t evslot)
{
// Reuse the next track booked for a given event slot.

   fTracksLock.Lock();

   Int_t itrack = fMaxPerEvent*evslot+fNtracks[evslot];
   GeantTrack *track = fTracks[itrack];
   track->Reset();
   track->particle = itrack;
   fEvents[evslot]->AddTrack();
   fNtracks[evslot]++;
//   printf("slot %d, track %d = %p\n", evslot, fNtracks[evslot], track);

   if (fNtracks[evslot]==fMaxPerEvent) {
      Fatal("AddTrack", "No room to add track");
   }

   fTracksLock.UnLock();

   return track;
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(GeantTrack *track)
{
// Mark track as stopped for tracking.
//   Printf("Stopping track from event %d in slot %d. Left %d tracks.", fEvents[track->evslot]->event, track->evslot, fEvents[track->evslot]->TracksLeft());
   if (track->IsAlive()) {
      fEvents[track->evslot]->StopTrack();
      fNtransported++;
   } else {
      Printf ("Stopping dead track.");
   }
   track->Kill();

	if (fEvents[track->evslot]->Transported()) {
      PerThread::reference TBBperThread = fTBBthreadData.local();
		TBBperThread.PushFinishedEvent (track->evslot);

//      Printf ("Finished event %d in slot %d\tfNcurPriorEvents=%d", track->event, track->evslot, fNcurPriorEvents);
      Printf ("Finished event %d in slot %d", track->event, track->evslot);

      fEventsStatus[fEvents[track->evslot]->event]++;    // In principle this can be moved into the PropTask

	}

}

//______________________________________________________________________________
Int_t GeantPropagator::ImportTracks(Int_t nevents, Double_t average, Int_t startevent, Int_t startslot)
{
// Import tracks from "somewhere". Here we just generate nevents.
   PerThread::reference TBBperThread = fTBBthreadData.local();

   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
   *TBBperThread.fMatrix = nav->GetCurrentMatrix();
   TBBperThread.fVolume = node->GetVolume();
   GeantVolumeBasket *basket = (GeantVolumeBasket*)TBBperThread.fVolume->GetField();
   TGeoBranchArray a;
   a.InitFromNavigator(gGeoManager->GetCurrentNavigator());

   const Double_t etamin = -3, etamax = 3;
   Int_t ntracks = 0;
   Int_t ntracksTotal = 0;

   // Species generated for the moment N, P, e, photon
   const Int_t kMaxPart=9;
   const Int_t pdgGen[9] =        {kPiPlus, kPiMinus, kProton, kProtonBar, kNeutron, kNeutronBar, kElectron, kPositron, kGamma};
   const Double_t pdgRelProb[9] = {   1.,       1.,      1.,        1.,       1.,          1.,        1.,        1.,     1.};
   const Species_t pdgSpec[9] =    {kHadron, kHadron, kHadron, kHadron, kHadron, kHadron, kLepton, kLepton, kLepton};
   static Double_t pdgProb[9] = {0.};
   Int_t pdgCount[9] = {0};

   static Bool_t init=kTRUE;
   if(init) {
      pdgProb[0]=pdgRelProb[0];
      for(Int_t i=1; i<kMaxPart; ++i) pdgProb[i]=pdgProb[i-1]+pdgRelProb[i];
      init=kFALSE;
   }
   Int_t event = startevent;
   for (Int_t slot=startslot; slot<startslot+nevents; slot++) {
      ntracks = TBBperThread.fRndm->Poisson(average);

      ntracksTotal += ntracks;

      Printf ("Importing event %d into slot %d with %d tracks.", event, slot, ntracks);

      if (!fEvents[slot]) fEvents[slot] = new GeantEvent();
      fEvents[slot]->slot = slot;
      fEvents[slot]->event = event;
      fEvents[slot]->Reset();
      fNtracks[slot] = 0;

      for (Int_t i=0; i<ntracks; i++) {
         GeantTrack *track = AddTrack(slot);
         *track->path = a;
         *track->nextpath = a;
         track->event = event;
         track->evslot = slot;
         Double_t prob=TBBperThread.fRndm->Uniform(0.,pdgProb[kMaxPart-1]);
         track->pdg=0;
         for(Int_t j=0; j<kMaxPart; ++j) {
            if(prob <= pdgProb[j]) {
               track->pdg = pdgGen[j];
               track->species = pdgSpec[j];
//            Printf("Generating a %s",TDatabasePDG::Instance()->GetParticle(track->pdg)->GetName());
               pdgCount[j]++;
               break;
            }
         }
         if(!track->pdg) Fatal("ImportTracks","No particle generated!");
         TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track->pdg);
         track->charge = part->Charge()/3.;
         track->mass   = part->Mass();
         track->xpos = fVertex[0];
         track->ypos = fVertex[1];
         track->zpos = fVertex[2];
         track->e = fKineTF1->GetRandom()+track->mass;
         Double_t p = TMath::Sqrt((track->e-track->mass)*(track->e+track->mass));
         Double_t eta = TBBperThread.fRndm->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*TBBperThread.fRndm->Rndm();
         track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
         track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
         track->pz = p*TMath::Cos(theta);
         track->frombdr = kFALSE;
         Int_t itrack = track->particle;
	
         TBBperThread.fCollection->AddTrack(itrack, basket);

         fNstart++;
      }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      event++;
      for (Int_t i=0; i<kMaxPart; i++) {
//         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
         pdgCount[i] = 0;
      }

   }

   fNadded += ntracksTotal;

//   Printf("Injecting %d events...", nevents);
   InjectCollection(TBBperThread.fCollection);
   return ntracksTotal;
}


//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance()
{
// Single instance of the propagator
   if (!fgInstance) fgInstance = new GeantPropagator();
   return fgInstance;
}

//______________________________________________________________________________
void GeantPropagator::Initialize()
{
// Initialization
   fMaxPerEvent = 5*fNaverage;
   fMaxTracks = fMaxPerEvent*fNevents;
	
// Initialize arrays here.
   if (!fKineTF1) {
      fKineTF1 = new TF1("fKineTF1","gaus",fEmin,fEmax);
      fKineTF1->SetParameters(1,3*fEmin,5);
   }

   if (!fProcesses) {
      fProcesses = new PhysicsProcess*[fNprocesses];
      fProcesses[0] = new ScatteringProcess("Scattering");
      fProcesses[1] = new ElossProcess("Eloss");
      fElossInd = 1;
      fProcesses[2] = new InteractionProcess("Interaction");
   }

   if(!fNtracks){
     fNtracks = new Int_t[fNevents];
     memset(fNtracks,0,fNevents*sizeof(Int_t));
   }
   if (!fTracks) {
      fTracks = new GeantTrack*[fMaxTracks];
      memset(fTracks, 0, fMaxTracks*sizeof(GeantTrack*));
      // Pre-book tracks (Thanks to olav.smorholm@cern.ch)
      fTracksStorage.reserve(fMaxTracks);
      for (Int_t i=0; i<fMaxTracks; i++) {
         fTracksStorage.push_back(GeantTrack(0));
         fTracks[i] = &(fTracksStorage.back());
      }
   }

   GeantThreadData td(fMaxPerBasket, 3);
   fTBBthreadData = PerThread(td);

   fDispThrPriority = 50;
   fEventsStatus = new tbb::atomic<Int_t>[fNtotal];
   for (Int_t evn=0; evn<fNtotal; evn++) fEventsStatus[evn] = 0;
   fPriorityRange[0] = -1;
   fPriorityRange[1] = -1;
   fPrioritize = kFALSE;

   niter = -1;

   fWMgr = WorkloadManager::Instance();
   // Add some empty baskets in the queue
   fWMgr->AddEmptyBaskets(1000);


	numOfPNtasks = new TH1S ("numOfPNtasks", "Number of normal PropTasks running; (1/5)second", 1500, 0, 1500);
	numOfPPtasks = new TH1S ("numOfPPtasks", "Number of priority PropTasks running; (1/5)second", 1500, 0, 1500);
	numOfDtasks = new TH1S ("numOfDtasks", "Number of DispTasks running; (1/5)second", 1500, 0, 1500);
	sizeOfFQ = new TH1S ("sizeOfFQ", "Size of feeder_queue; (1/5)second", 1500, 0, 1500);
	sizeOfPFQ = new TH1S ("sizeOfPFQ", "Size of priority_feeder_queue; (1/5)second", 1500, 0, 1500);
	sizeOfCQ = new TH1S ("sizeOfCQ", "Size of collector_queue; (1/5)second", 1500, 0, 1500);
	numOfCollsPerTask = new TH1S ("numOfCollsPerTask", "Number of CollDispTasks; Number of collections poped", 60, 0, 60);
	numOfTracksInBasket = new TH1S ("numOfTracksInBasket", "Number of normal baskets; tracks per basket", fNperBasket*2, 0, fNperBasket*2);
	numOfTracksInPriorBasket = new TH1S ("numOfTracksInPriorBasket", "Number of priority baskets; tracks per basket", fNperBasket*2, 0, fNperBasket*2);
	numOfTracksInColl = new TH1S ("numOfTracksInColl", "Number of collection; tracks per collection", 2 * fNaverage, 0, 2 * fNaverage);
   numOfTracksTransportedInTime = new TH1I ("numOfTracksTransportedInTime", "Tracks transported; (1/5)second", 1500, 0, 1500);

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

//______________________________________________________________________________
void GeantPropagator::InjectCollection(GeantTrackCollection* inColl)
{
// Inject collector handled by a single thread in the collector queue.

   GeantTrackCollection *newcoll;

   if (fWMgr->tbb_collector_empty_queue.try_pop(newcoll)) { }
   else newcoll = new GeantTrackCollection(100);

   *newcoll = *inColl;
   fWMgr->tbb_collector_queue.push(newcoll);
   inColl->Clear();
}

//______________________________________________________________________________
Bool_t GeantPropagator::LoadGeometry(const char *filename)
{
// Load the detector geometry from file.
   if (gGeoManager) return kTRUE;
   TGeoManager *geom = TGeoManager::Import(filename);
   if (geom) {
      // Create the basket array
      Int_t nvols = gGeoManager->GetListOfVolumes()->GetEntries();
      fWMgr->CreateBaskets(2*nvols);
      return kTRUE;
   }
   ::Error("LoadGeometry","Cannot load geometry from file %s", filename);
   return kFALSE;
}

//______________________________________________________________________________
void GeantPropagator::PhysicsSelect(Int_t ntracks, Int_t *trackin)
{
// Generate all physics steps for the tracks in trackin.
// Vectorized, except the unavoidable Sort()

   PerThread::reference TBBperThread = fTBBthreadData.local();

   static const Double_t maxlen = TMath::Limits<double>::Max();
   Double_t pstep;
   Int_t ipart, iproc;
   GeantTrack *track;
   Double_t *procStep;
   // Fill interaction lengths for all processes and all particles
   for (iproc=0; iproc<fNprocesses; iproc++) {
      procStep = TBBperThread.GetProcStep(iproc);
      fProcesses[iproc]->ComputeIntLen(TBBperThread.fVolume, ntracks, trackin, procStep);
   }
   // Loop tracks and select process
   for (Int_t i=0; i<ntracks; i++) {
      ipart = trackin[i];
      track = fTracks[ipart];
      track->step = maxlen;
      track->process = -1;
      for (iproc=0; iproc<fNprocesses; iproc++) {
         procStep = TBBperThread.GetProcStep(iproc);
         pstep = procStep[i];
         if (pstep < track->step) {
            track->step = pstep;
            track->process = iproc;
         }
      }
      if (fUseDebug && (fDebugTrk==ipart || fDebugTrk<0)) {
         Printf("   PhysicsSelect: track #%d - process=%d pstep=%g",ipart,track->process,track->step);
         if (track->step>1.E200) {
            Printf("xxx");
         }
      }
   }
}

//______________________________________________________________________________
void GeantPropagator::PrintParticles(Int_t *trackin, Int_t ntracks)
{
// Print the detailed particles list.
   Printf("================ THREAD : particles list");
   for (Int_t i=0; i<ntracks; i++) {
      fTracks[trackin[i]]->Print();
   }
}

//______________________________________________________________________________
void* GeantPropagator::GlobalObserver (void* arg)
{
	Printf("======== Observer started ========");
	Int_t recvSig;

	GeantPropagator* prop = GeantPropagator::Instance ();
	Int_t i=0;

	while (1)
	{
		if (prop->observerSigsQueue.try_pop(recvSig))
		{
			break;
		}
		gSystem->Sleep (200);
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
void GeantPropagator::PropagatorGeom(const char *geomfile, Bool_t graphics, Bool_t single)
{
// Propagate fNevents in the volume containing the vertex.
// Simulate 2 physics processes up to exiting the current volume.
   static Bool_t called=kFALSE;
   fUseGraphics = graphics;
   fSingleTrack = single;
   Initialize();
   if (called) {
      Printf("Sorry, you can call this only once per session.");
      return;
   }
   called = kTRUE;

   // Initialize geometry and current volume
   if (!LoadGeometry(geomfile)) return;
   if (fSingleTrack) Printf("==== Executing in single track loop mode using %d threads ====", fNthreads);
   else              Printf("==== Executing in vectorized mode using %d threads ====",fNthreads);
   if (fFillTree)    Printf("  I/O enabled - disable if comparing single track loop with vectorized modes");
   else              Printf("  I/O disabled");
   if (fUsePhysics)  Printf("  Physics ON with %d processes", fNprocesses);
   else              Printf("  Physics OFF");
   // Create the basket array

   // Import the input events. This will start also populating the main queue
   if (!fEvents) {
      fEvents = new GeantEvent*[fNevents];
      memset(fEvents, 0, fNevents*sizeof(GeantEvent*));
   }

   Int_t NinjectedTracks = ImportTracks(fNevents, fNaverage, 0, 0);
   fNimportedEvents += fNevents;

   // Initialize tree
   fOutput = new GeantOutput();
   fOutput->Init(fMaxTracks);
   if (fFillTree) {
      fOutFile = new TFile("output.root", "RECREATE");
      fOutTree = new TTree("TK","Transport track data");
      fOutTree->Branch("gen", &fOutput);
   }

   // Loop baskets and transport particles until there is nothing to transport anymore
   gGeoManager->SetMaxThreads(fNthreads+1);
   tbb::task_scheduler_init init(fNthreads+1);

   fTimer = new TStopwatch();

/*------------------------------------------------------------------------------------------------------*/
// Just histograms filling
///*
   if (fUseGraphics) {
      hnb = new TH1F("hnb","number of baskets in the transport queue; iteration#",200000,0,200000);
      hnb->SetMinimum(0.);
      hnb->SetMaximum(NinjectedTracks / fNperBasket + 200.);
      hnb->SetFillColor(kRed);
      htracks = new TH1F("htracks","number of tracks/basket; iteration#",200000,0,200000);
      htracks->SetFillColor(kGreen);
      numOfTracksTransportedInIter = new TH1I("numOfTracksTransportedInIter","number of tracks transported; iteration#",200000,0,200000);
      numOfTracksTransportedInIter->SetFillColor(kBlue);
   }
//*/
/*------------------------------------------------------------------------------------------------------*/

	TThread* observerThread = new TThread (GeantPropagator::GlobalObserver);
	observerThread->Run();

   fTimer->Start();

   CollDispTask& topTask = *new (task::allocate_root()) CollDispTask(1);
   task::spawn_root_and_wait (topTask);

   // Here we check if really all the events were transported
   Int_t first_not_transported = fNtotal;
   for (Int_t evn=0; evn<fNtotal; evn++) {
      if (fEventsStatus[evn] == 0) {
         first_not_transported = evn;
         break;
      }
   }
   if (first_not_transported < fNtotal) {
      Printf ("Not all events transported. First not transported: %d", first_not_transported);
      Printf ("Feeder queue size: %d\t\tCollector queue size: %d",
               fWMgr->tbb_feeder_queue.size(), fWMgr->tbb_collector_queue.size());
      Printf ("Collections waiting: %d\t\tTracks waiting: %d", fCollsWaiting, fTracksWaiting);

      fCollsWaiting = 0;
      fTracksWaiting = 0;

      fGarbageCollMode = kTRUE;
      CollDispTask& garbTask = *new (task::allocate_root()) CollDispTask(fWMgr->tbb_collector_queue.size());
      task::spawn_root_and_wait (garbTask);
   }

   first_not_transported = fNtotal;
   for (Int_t evn=0; evn<fNtotal; evn++) {
      if (fEventsStatus[evn] == 0) {
         first_not_transported = evn;
         break;
      }
   }
   if (first_not_transported < fNtotal) Printf ("Not all events transported. First not transported: %d", first_not_transported);
   else Printf ("ALL EVENTS FINISHED");

   fTimer->Stop();

   observerSigsQueue.push(1);

/*------------------------------------------------------------------------------------------------------*/
// Just histograms filling
///*
   TPad *pad1=0, *pad2=0, *pad3=0;
   TCanvas *c1 = 0;
   if (fUseGraphics) {
      c1 = new TCanvas("c2","c2",1200,1200);
      c1->Divide(1,3);

      pad1 = (TPad*)c1->cd(1);
      hnb->SetStats(kFALSE);
      hnb->GetXaxis()->SetRangeUser(0,niter/10);
      hnb->Draw();

      pad2 = (TPad*)c1->cd(2);
      htracks->SetStats(kFALSE);
      htracks->GetXaxis()->SetRangeUser(0,niter/10);
      htracks->Draw();

      pad3 = (TPad*)c1->cd(3);
      numOfTracksTransportedInIter->SetStats(kFALSE);
      numOfTracksTransportedInIter->GetXaxis()->SetRangeUser(0,niter/10);
      numOfTracksTransportedInIter->Draw();
   }
//*/
/*------------------------------------------------------------------------------------------------------*/

	TFile* diagFile = new TFile("diagnostics.root", "UPDATE");

   hnb->Write();
   htracks->Write();
   numOfTracksTransportedInIter->Write();

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

/*------------------------------------------------------------------------------------------------------*/

   Double_t rtime = fTimer->RealTime();
   Double_t ctime = fTimer->CpuTime();
   if (fFillTree) fOutTree->AutoSave();
   delete fOutFile;

   FILE *fp = fopen("rez.txt","a");
   fprintf(fp, "threads %d\nevtotal %d\nevbuf %d\ntr/ev %f\ntr/bask %d\n",
               fNthreads, fNtotal, fNevents, fNaverage, fNperBasket);
   fprintf(fp, "added %lld\ntransported %lld\nsafesteps %lld\nsnextsteps %lld\nrealtime %g\ncputime %g\n",
               fNadded, fNtransported, fNsafeSteps, fNsnextSteps, rtime, ctime);
   fclose(fp);

   fOutFile = 0;
   fOutTree = 0;
}

//______________________________________________________________________________
void GeantPropagator::SelectTracksForProcess(Int_t iproc, Int_t ntotransport, Int_t *particles, Int_t &ntodo, Int_t *parttodo)
{
// Add to output array all particles that will do the process iproc.
   for (Int_t itr=0; itr<ntotransport; itr++)
      if (fTracks[particles[itr]]->process == iproc) parttodo[ntodo++] = particles[itr];
}

