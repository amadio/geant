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
#include "GeantOutput.h"
#include "sync_objects.h"
#include "PhysicsProcess.h"
#include "GeantVolumeBasket.h"
#include "WorkloadManager.h"
#include "GeantThreadData.h"

GeantPropagator *gPropagator = 0;
   
ClassImp(GeantPropagator)

GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
                :TObject(),
                 fNthreads(1),
                 fNevents(100),
                 fNtotal(1000),
                 fNtracks(0),
                 fNtransported(0),
                 fNsafeSteps(0),
                 fNsnextSteps(0),
                 fNprocesses(3),
                 fElossInd(0),
                 fNstart(0),
                 fMaxTracks(0),  //change
                 fMaxThreads(100),
                 fNminThreshold(10),
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
                 fTransportOngoing(kFALSE),
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
                 fTracksPerBasket(0),
                 fWaiting(0),
                 fThreadData(0)
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
   
   if (fThreadData) {
      for (i=0; i<fNthreads; i++) delete fThreadData[i];
      delete [] fThreadData;
   }   
   delete [] fWaiting;
   delete fOutput;
   delete fOutFile;
   delete fTimer;
   delete fWMgr;
}

//______________________________________________________________________________
Int_t GeantPropagator::AddTrack(GeantTrack *track)
{
// Add a new track in the system.
   TThread::Lock();
   Int_t slot = track->evslot;
   Int_t itrack = fMaxPerEvent*slot+fNtracks[slot];
   track->particle = itrack;
   fEvents[slot]->AddTrack();
//   Int_t tid = TGeoManager::ThreadId();
   fTracks[itrack] = track;
   fNtracks[slot]++;
   fNtransported++;
   if (fNtracks[slot]==fMaxPerEvent) {
      Fatal("AddTrack", "No room to add track");
   }
   TThread::UnLock();
   return itrack;
}

//______________________________________________________________________________
GeantTrack_v *GeantPropagator::AddTrack(Int_t evslot, const TGeoBranchArray &start)
{
// Add a track for the event stored at evslot. The tracks starts in the specified
// location.
   Int_t itrack = fMaxPerEvent*evslot+fNtracks[evslot];
   GeantTrack *track = fTracks[itrack];
   track->Reset();
   track->particle = itrack;
   fEvents[evslot]->AddTrack();
   fTracksLock.Lock();
// critical section
   fNtracks[evslot]++;
   fNtransported++;
//   printf("slot %d, track %d = %p\n", evslot, fNtracks[evslot], track);
// end critical section   
   fTracksLock.UnLock();
   if (fNtracks[evslot]==fMaxPerEvent) {
      Fatal("AddTrack", "No room to add track");
   }
   return track;
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(GeantTrack *track)
{
// Mark track as stopped for tracking.
//   Printf("Stopping track %d", track->particle);
   if (track->IsAlive()) fEvents[track->evslot]->StopTrack();
   track->Kill();
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack()
{
// Returns a temporary track support for the physics processes, unique per
// thread which can be used to add tracks produced by physics processes.
   Int_t tid = TGeoManager::ThreadId();
   return fThreadData[tid]->fTrack;
}
   
//______________________________________________________________________________
GeantVolumeBasket *GeantPropagator::ImportTracks(Int_t nevents, Double_t average, Int_t startevent, Int_t startslot)
{
// Import tracks from "somewhere". Here we just generate nevents.
   Int_t tid = TGeoManager::ThreadId();
   GeantThreadData *td = fThreadData[tid];
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();
   TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
   *td->fMatrix = nav->GetCurrentMatrix();
   td->fVolume = node->GetVolume();
   GeantVolumeBasket *basket = (GeantVolumeBasket*)td->fVolume->GetField();
   TGeoBranchArray a;
   a.InitFromNavigator(gGeoManager->GetCurrentNavigator());
   
   const Double_t etamin = -3, etamax = 3;
   Int_t ntracks = 0;
   
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
      ntracks = td->fRndm->Poisson(average);
      if (!fEvents[slot]) fEvents[slot] = new GeantEvent();
      fEvents[slot]->slot = slot;
      fEvents[slot]->event = event;
      fEvents[slot]->Reset();
      fNtracks[slot] = 0;
      
      for (Int_t i=0; i<ntracks; i++) {
//         TGeoBranchArray *a = new TGeoBranchArray();
//         a->InitFromNavigator(gGeoManager->GetCurrentNavigator());
         GeantTrack *track = AddTrack(slot);
         *track->path = a;
         *track->nextpath = a;
         track->event = event;
         track->evslot = slot;
         Double_t prob=td->fRndm->Uniform(0.,pdgProb[kMaxPart-1]);
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
         Double_t eta = td->fRndm->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*td->fRndm->Rndm();
         track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
         track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
         track->pz = p*TMath::Cos(theta);
         track->frombdr = kFALSE;
         Int_t itrack = track->particle;
	 
         fCollections[tid]->AddTrack(itrack, basket);
//         gPropagator->fCollections[tid]->AddTrack(itrack, basket);
    //     basket->AddTrack(fNstart);
         fNstart++;
      }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      event++;
      for (Int_t i=0; i<kMaxPart; i++) {
//         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
         pdgCount[i] = 0;
      }   
   }
//   Printf("Injecting %d events...", nevents);
   InjectCollection(tid);      
   return basket;
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
   gPropagator = GeantPropagator::Instance();
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
   
   if (!fTracksPerBasket) {
      fTracksPerBasket = new Int_t[fNthreads];    
      for (Int_t i=0; i<fNthreads; i++) fTracksPerBasket[i] = 0;
   }   
   if (!fCollections) {
      fCollections = new GeantTrackCollection*[fNthreads+1];    
      for (Int_t i=0; i<fNthreads+1; i++) fCollections[i] = new GeantTrackCollection(100);
   }
   if (!fWaiting) {
      fWaiting = new UInt_t[fNthreads+1];
      memset(fWaiting, 0, (fNthreads+1)*sizeof(UInt_t));
   }  
   if (!fThreadData) {
      fThreadData = new GeantThreadData*[fNthreads];
      for (Int_t i=0; i<fNthreads+1; i++) fThreadData[i] = new GeantThreadData(fMaxPerBasket, 3);
   } 
   fWMgr = WorkloadManager::Instance(fNthreads);
   // Add some empty baskets in the queue
   fWMgr->AddEmptyBaskets(1000);
}

//______________________________________________________________________________
UInt_t GeantPropagator::GetNwaiting() const
{
// Returns number of waiting threads. Must be called by a single scheduler.
   UInt_t nwaiting = 0;
   for (Int_t tid=0; tid<fNthreads; tid++) nwaiting += fWaiting[tid];
   return nwaiting;
}

//______________________________________________________________________________
GeantBasket *GeantPropagator::InjectBasket(GeantBasket *basket)
{
// Inject basket in the work queue and return a fresh one.
   if (!basket->GetNtracks()) return basket;
   GeantBasket *newbasket;
   if (fWMgr->EmptyQueue()->empty()) newbasket = new GeantBasket(fNperBasket);
   else newbasket = (GeantBasket*)fWMgr->EmptyQueue()->wait_and_pop();
   fWMgr->FeederQueue()->push(basket);
   return newbasket;
}

//______________________________________________________________________________
void GeantPropagator::InjectCollection(Int_t tid)
{
// Inject collector handled by a single thread in the collector queue.
//   if (!fCollections[tid]->GetNtracks()) return;
   GeantTrackCollection *newcoll;
   if (fWMgr->CollectorEmptyQueue()->empty()) newcoll = new GeantTrackCollection(100);
   else newcoll = (GeantTrackCollection*)fWMgr->CollectorEmptyQueue()->wait_and_pop();
   GeantTrackCollection *toinject = fCollections[tid];
   fCollections[tid] = newcoll;
   fWMgr->CollectorQueue()->push(toinject);
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
void GeantPropagator::PhysicsSelect(Int_t ntracks, Int_t *trackin, Int_t tid)
{
// Generate all physics steps for the tracks in trackin.
// Vectorized, except the unavoidable Sort()
   static const Double_t maxlen = TMath::Limits<double>::Max();   
   Double_t pstep;
   Int_t ipart, iproc;
   GeantTrack *track;
   GeantThreadData *td = fThreadData[tid];
   Double_t *procStep;
   // Fill interaction lengths for all processes and all particles
   for (iproc=0; iproc<fNprocesses; iproc++) {
      procStep = td->GetProcStep(iproc);
      fProcesses[iproc]->ComputeIntLen(td->fVolume, ntracks, trackin, procStep, tid);
   }   
   // Loop tracks and select process
   for (Int_t i=0; i<ntracks; i++) {
      ipart = trackin[i];
      track = fTracks[ipart];
      track->step = maxlen;
      track->process = -1;
      for (iproc=0; iproc<fNprocesses; iproc++) {
         procStep = td->GetProcStep(iproc);
         pstep = procStep[i];
         if (pstep < track->step) {
            track->step = pstep;
            track->process = iproc;
         }
      }
      if (fUseDebug && (fDebugTrk==ipart || fDebugTrk<0)) {
         Printf("   (%d) PhysicsSelect: track #%d - process=%d pstep=%g",tid,ipart,track->process,track->step);
         if (track->step>1.E200) {
            Printf("xxx");
         }   
      }
   }      
}

//______________________________________________________________________________
void GeantPropagator::PrintParticles(Int_t *trackin, Int_t ntracks, Int_t tid)
{
// Print the detailed particles list.
   Printf("================ THREAD %d: particles list", tid);
   for (Int_t i=0; i<ntracks; i++) {
      fTracks[trackin[i]]->Print();
   }
}
      
//______________________________________________________________________________
void GeantPropagator::PropagatorGeom(const char *geomfile, Int_t nthreads, Bool_t graphics, Bool_t single)
{
// Propagate fNevents in the volume containing the vertex. 
// Simulate 2 physics processes up to exiting the current volume.
   static Bool_t called=kFALSE;
   fUseGraphics = graphics;
   fNthreads = nthreads;
   fSingleTrack = single;
   Initialize();
   if (called) {
      Printf("Sorry, you can call this only once per session.");
      return;
   }
   called = kTRUE;
   
//   Int_t itrack;

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
   
   ImportTracks(fNevents, fNaverage, 0, 0);

   // Initialize tree
   fOutput = new GeantOutput();
   fOutput->Init(fMaxTracks);
   if (fFillTree) {
      fOutFile = new TFile("output.root", "RECREATE");
      fOutTree = new TTree("TK","Transport track data");
      fOutTree->Branch("gen", &fOutput);
   }
      
   // Loop baskets and transport particles until there is nothing to transport anymore
   fTransportOngoing = kTRUE;
   gGeoManager->SetMaxThreads(nthreads);
   fTimer = new TStopwatch();
   fWMgr->StartThreads();
   fTimer->Start();
   fWMgr->WaitWorkers();
   fTimer->Stop();
   Double_t rtime = fTimer->RealTime();
   Double_t ctime = fTimer->CpuTime();
   if (fFillTree) fOutTree->AutoSave();
   delete fOutFile;
//   fTimer->Print();
   Double_t speedup = ctime/rtime;
   Double_t efficiency = speedup/nthreads;
//   fWMgr->Print();
   fWMgr->JoinThreads();
   const char *geomname=geomfile;
   if(strstr(geomfile,"http://root.cern.ch/files/")) geomname=geomfile+strlen("http://root.cern.ch/files/");
   Printf("=== Transported: %lld,  safety steps: %lld,  snext steps: %lld, RT=%gs, CP=%gs", fNtransported, fNsafeSteps, fNsnextSteps,rtime,ctime);
   Printf("   nthreads=%d + 1 garbage collector speed-up=%f  efficiency=%f", nthreads, speedup, efficiency);
   gSystem->mkdir("results");
   FILE *fp = fopen(Form("results/%s_%d.dat",geomname,single),"w");
   fprintf(fp,"%d %lld %lld %g %g",single, fNsafeSteps, fNsnextSteps,rtime,ctime);
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

//______________________________________________________________________________
Bool_t DrawData(Int_t color = kRed)
{
// Draw a given generation of track points.
   static Long64_t ientry = 0;
   Long64_t nentries = gPropagator->fOutTree->GetEntries();
   if (ientry==nentries) return kFALSE;   
   TPolyMarker3D *pmgen = new TPolyMarker3D();
   pmgen->SetMarkerColor(color);
   Int_t nread = 0;
   Int_t ibgen = 0;
   while (ientry<nentries) {
      gPropagator->fOutTree->GetEntry(ientry++);
      if (!nread++) {
         ibgen = gPropagator->fOutput->fBasketGeneration;
      }   
      if (gPropagator->fOutput->fBasketGeneration > ibgen) {
         ientry--;
         break;
      }      
      for (Int_t itrack=0; itrack<gPropagator->fOutput->fNtracks; itrack++) 
         pmgen->SetNextPoint(gPropagator->fOutput->fX[itrack], gPropagator->fOutput->fY[itrack],gPropagator->fOutput->fZ[itrack]);
   }
   Printf("basket generation #%d\n", ibgen);
   pmgen->Draw("SAME");
   if (ientry==nentries) return kFALSE;
   return kTRUE;
}

//______________________________________________________________________________
void DrawNextBasket()
{
// Draw next basket
   Bool_t drawn = kTRUE;
   if (!gPropagator->fOutFile) gPropagator->fOutFile = new TFile("output.root");
   if (!gPropagator->fOutTree) {
      gPropagator->fOutTree = (TTree*)gPropagator->fOutFile->Get("TK");
      gPropagator->fOutput = new GeantOutput();
      gPropagator->fOutTree->SetBranchAddress("gen", &gPropagator->fOutput);
   }   
   TProcessEventTimer *timer = new TProcessEventTimer(1);
   gROOT->SetInterrupt(kFALSE);
   
   while (drawn) {
      if (gROOT->IsInterrupted()) break;
      if (timer->ProcessEvents()) continue;
      drawn = DrawData(Int_t(8*gRandom->Rndm())+1);
      gPad->Modified();
      gPad->Update();
      if (!drawn) {
         Printf("That was the last basket...\n");
      }
   }
   gROOT->SetInterrupt(kTRUE);      
}   

//______________________________________________________________________________
void Stop()
{
// Stop the process timer loop.
   gROOT->SetInterrupt(kTRUE);
}   

//______________________________________________________________________________
void Menu(const char *file="geometry.root")
{
// Start a TControl bar menu
   GeantPropagator *propagator = GeantPropagator::Instance();
   if (!gGeoManager) propagator->LoadGeometry(file);
   gGeoManager->SetVisLevel(1);
   gGeoManager->GetTopVolume()->Draw("ogl");
   TGLSAViewer *viewer = (TGLSAViewer *)gPad->GetViewer3D();
   viewer->SetResetCamerasOnUpdate(kFALSE);
   TControlBar *bar = new TControlBar("vertical", "Propagator",100,10);
   bar->AddButton("Run","PropagatorGeom()", "You can run this only once"); 
   bar->AddButton("ShowTracks", "DrawNextBasket()", "Draw next generation of baskets");
   bar->AddButton("Stop", "Stop()", "Stop drawing.");
   bar->Show();
   gROOT->SaveContext();
   Printf("=== Maybe press 'w' for wireframe mode ===\n");
}
