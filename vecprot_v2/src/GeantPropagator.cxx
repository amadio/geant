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
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TGenPhaseSpace.h"
#include "GeantTrack.h"
#include "GeantOutput.h"
#include "sync_objects.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"
#include "GeantBasket.h"
#include "GeantThreadData.h"
#include "GeantVApplication.h"
#include "GeantFactoryStore.h"
#include "GeantEvent.h"
#include "GeantScheduler.h"

GeantPropagator *gPropagator = 0;
   
ClassImp(GeantPropagator)

GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
                :TObject(),
                 fNthreads(1),
                 fNevents(100),
                 fNtotal(1000),
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
                 fMaxPerBasket(10000),
                 fMaxPerEvent(0),
                 fNaverage(0.),
                 fVertex(),
                 fEmin(1.E-4), // 100 KeV
                 fEmax(10),  // 10 Gev
                 fBmag(1.),
                 fUsePhysics(kTRUE),
                 fUseDebug(kFALSE),
                 fUseGraphics(kFALSE),
                 fTransportOngoing(kFALSE),
                 fSingleTrack(kFALSE),
                 fFillTree(kFALSE),
                 fTracksLock(),
                 fWMgr(0),
                 fApplication(0),
                 fOutput(0),
                 fOutTree(0),
                 fOutFile(0),
                 fTimer(0),
//                 fProcesses(0),
                 fProcess(0),
                 fStoredTracks(0),
                 fNtracks(0),
                 fEvents(0),
                 fWaiting(0),
                 fThreadData(0)
{
// Constructor
   fVertex[0] = -8.;
   fVertex[1] = fVertex[2] = 0.;
//   for (Int_t i=0; i<3; i++) fVertex[i] = gRandom->Gaus(0.,0.2);
   fgInstance = this;
}

//______________________________________________________________________________
GeantPropagator::~GeantPropagator()
{
// Destructor
   Int_t i;
   delete fProcess;
//   if (fProcesses) {
//     for (i=0; i<fNprocesses; i++) delete fProcesses[i];
//     delete [] fProcesses;
//   }  

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
   delete fApplication;
}

//______________________________________________________________________________
Int_t GeantPropagator::AddTrack(GeantTrack &track)
{
// Add a new track in the system. returns track number within the event.
   Int_t slot = track.fEvslot;
   track.fParticle = fEvents[slot]->AddTrack();
//   fNtracks[slot]++;
//   fNtransported++;
   return track.fParticle;
}

//______________________________________________________________________________
Int_t GeantPropagator::DispatchTrack(const GeantTrack &track)
{
// Dispatch a registered track produced by the generator.
   return fWMgr->GetScheduler()->AddTrack(track);
}   
   
//______________________________________________________________________________
void GeantPropagator::StopTrack(GeantTrack *track)
{
// Mark track as stopped for tracking.
//   Printf("Stopping track %d", track->particle);
   if (track->IsAlive()) fEvents[track->fEvslot]->StopTrack();
   track->Kill();
}

//______________________________________________________________________________
void GeantPropagator::StopTrack(const GeantTrack_v &tracks, Int_t itr)
{
// Mark track as stopped for tracking.
//   Printf("Stopping track %d", track->particle);
   fEvents[tracks.fEvslotV[itr]]->StopTrack();
}

//______________________________________________________________________________
GeantTrack &GeantPropagator::GetTempTrack(Int_t tid)
{
// Returns a temporary track support for the physics processes, unique per
// thread which can be used to add tracks produced by physics processes.
   if (tid<0) tid = TGeoManager::ThreadId();
   GeantTrack &track = fThreadData[tid]->fTrack;
   track.Clear();
   track.fPath = new TGeoBranchArray();
   track.fNextpath = new TGeoBranchArray();
   return track;
}
   
//______________________________________________________________________________
Int_t GeantPropagator::ImportTracks(Int_t nevents, Double_t average, Int_t startevent, Int_t startslot)
{
// Import tracks from "somewhere". Here we just generate nevents.
   static TGeoBranchArray *a = 0;
   Int_t tid = TGeoManager::ThreadId();
   GeantThreadData *td = fThreadData[tid];
   TGeoVolume *vol = 0;
   if (!a) {
      a = new TGeoBranchArray();
      TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
      if (!nav) nav = gGeoManager->AddNavigator();
      TGeoNode *node = nav->FindNode(fVertex[0], fVertex[1], fVertex[2]);
      *td->fMatrix = nav->GetCurrentMatrix();
      vol = node->GetVolume();
      td->fVolume = vol;
      a->InitFromNavigator(nav);
   } else {
      TGeoNode *node = a->GetCurrentNode();
      *td->fMatrix = a->GetMatrix();
      vol = node->GetVolume();
      td->fVolume = vol;
   }     
   GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr*>(vol->GetFWExtension());
   Int_t threshold = nevents*average/(2*fNthreads);
   threshold -= threshold%4;
   if (threshold<4) threshold = 4;
   if (threshold>256) threshold = 256;
   basket_mgr->SetThreshold(threshold);
   
   const Double_t etamin = -3, etamax = 3;
   Int_t ntracks = 0;
   Int_t ntotal = 0;
   Int_t ndispatched = 0;
   
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
      ntotal += ntracks;
      if (!fEvents[slot]) fEvents[slot] = new GeantEvent();
      fEvents[slot]->SetSlot(slot);
      fEvents[slot]->SetEvent(event);
      fEvents[slot]->Reset();
      
      for (Int_t i=0; i<ntracks; i++) {
         GeantTrack &track = GetTempTrack(tid);
         track.SetPath(a);
         track.SetNextPath(a);
         track.SetEvent(event);
         track.SetEvslot(slot);
         Double_t prob=td->fRndm->Uniform(0.,pdgProb[kMaxPart-1]);
//         track.SetPDG(kMuonPlus); // G5code=27
//         track.SetG5code(27);
         track.SetPDG(kElectron); // G5code=23
         track.SetG5code(23); // just a hack -> will change with new physics list
/*
         for(Int_t j=0; j<kMaxPart; ++j) {
            if(prob <= pdgProb[j]) {
               track.SetPDG(pdgGen[j]);
               track.SetSpecies(pdgSpec[j]);
//            Printf("Generating a %s",TDatabasePDG::Instance()->GetParticle(track->pdg)->GetName());
               pdgCount[j]++;
               break;
            }
         }   
*/
         if(!track.fPDG) Fatal("ImportTracks","No particle generated!");
         TParticlePDG *part = TDatabasePDG::Instance()->GetParticle(track.fPDG);
         track.SetCharge(part->Charge()/3.);
         track.SetMass(part->Mass());
         track.fXpos = fVertex[0];
         track.fYpos = fVertex[1];
         track.fZpos = fVertex[2];
         track.fE = 0.03 /*30MeV*/ +part->Mass();  //e-
//         track.fE = 0.3 /*300MeV*/ +part->Mass();  //mu+
         Double_t p = TMath::Sqrt((track.E()-track.Mass())*(track.E()+track.Mass()));
         track.SetP(p);
/*
         Double_t eta = td->fRndm->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*td->fRndm->Rndm();
         track.fXdir = TMath::Sin(theta)*TMath::Cos(phi);
         track.fYdir = TMath::Sin(theta)*TMath::Sin(phi);
         track.fZdir = TMath::Cos(theta);
*/
         track.fXdir = 1.;
         track.fYdir = 0.;
         track.fZdir = 0.;
         track.fFrombdr = kFALSE;
         track.fStatus = kAlive;
         
         AddTrack(track);
         ndispatched += DispatchTrack(track);
      }
//      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      event++;
      for (Int_t i=0; i<kMaxPart; i++) {
//         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
         pdgCount[i] = 0;
      }   
   }
   Printf("Imported %d tracks from events %d to %d. Dispatched %d baskets.",
           ntotal, startevent, startevent+nevents-1, ndispatched);
   return ndispatched;
}

//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance(Int_t ntotal, Int_t nbuffered)
{
// Single instance of the propagator
   if (!fgInstance) fgInstance = new GeantPropagator();
   if (ntotal) fgInstance->fNtotal = ntotal;
   if (nbuffered) {
      fgInstance->fNevents = nbuffered;
      GeantFactoryStore::Instance(nbuffered);
   }   
   return fgInstance;
}

/*
//______________________________________________________________________________
GeantPropagator *GeantPropagator::Instance()
{
// Single instance of the propagator
   if (!fgInstance) fgInstance = new GeantPropagator();
   return fgInstance;
}
*/
   
//______________________________________________________________________________
void GeantPropagator::Initialize()
{
// Initialization
   fMaxPerEvent = 5*fNaverage;  
   fMaxTracks = fMaxPerEvent*fNevents;
	
// Initialize arrays here.
   gPropagator = GeantPropagator::Instance();
   
   if (!fProcess) {
      Fatal("Initialize", "The physics process has to be initilaized before this");
      return;
   }
   // Initialize the process(es)
   fProcess->Initialize();
      
//   if (!fProcesses) {
//      fProcesses = new PhysicsProcess*[fNprocesses];
//      fProcesses[0] = new ScatteringProcess("Scattering");
//      fProcesses[1] = new ElossProcess("Eloss");
//      fElossInd = 1;
//      fProcesses[2] = new InteractionProcess("Interaction");
//   }

   if(!fNtracks){
     fNtracks = new Int_t[fNevents];
     memset(fNtracks,0,fNevents*sizeof(Int_t));
   }   
   
   if (!fWaiting) {
      fWaiting = new UInt_t[fNthreads+1];
      memset(fWaiting, 0, (fNthreads+1)*sizeof(UInt_t));
   }  
   if (!fThreadData) {
      fThreadData = new GeantThreadData*[fNthreads+1];
      for (Int_t i=0; i<fNthreads+1; i++) fThreadData[i] = new GeantThreadData(fMaxPerBasket, 3);
   } 
   fWMgr = WorkloadManager::Instance(fNthreads);
   // Add some empty baskets in the queue
//   fWMgr->AddEmptyBaskets(1000);
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
Bool_t GeantPropagator::LoadGeometry(const char *filename)
{
// Load the detector geometry from file.
   if (gGeoManager) return kTRUE;
   TGeoManager *geom = (gGeoManager)? gGeoManager : TGeoManager::Import(filename);
   if (geom) return kTRUE;
   ::Error("LoadGeometry","Cannot load geometry from file %s", filename);
   return kFALSE;
}

//______________________________________________________________________________
void GeantPropagator::PhysicsSelect(Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
{
// Generate all physics steps for the tracks in trackin.
   GeantThreadData *td = fThreadData[tid];
   // Reset the current step length to 0
   for (Int_t i=0; i<ntracks; ++i) {
      tracks.fStepV[i] = 0.;
      tracks.fEdepV[i] = 0.;
   }   
   fProcess->ComputeIntLen(td->fVolume->GetMaterial(), ntracks, tracks, 0, tid);

/*
   static const Double_t maxlen = TMath::Limits<double>::Max();   
   Double_t pstep;
   Int_t iproc;
   GeantThreadData *td = fThreadData[tid];
   Double_t *procStep;
   // Fill interaction lengths for all processes and all particles
   for (iproc=0; iproc<fNprocesses; iproc++) {
      if (fProcesses[iproc]->IsType(PhysicsProcess::kContinuous)) continue;
      procStep = td->GetProcStep(iproc);
      fProcesses[iproc]->ComputeIntLen(td->fVolume->GetMaterial(), ntracks, tracks, procStep, tid);
   }
   // Loop tracks and select process
   for (Int_t i=0; i<ntracks; i++) {
      // reset total step
      tracks.fStepV[i] = 0.;
      tracks.fPstepV[i] = maxlen;
      tracks.fProcessV[i] = -1;
      for (iproc=0; iproc<fNprocesses; iproc++) {
         if (fProcesses[iproc]->IsType(PhysicsProcess::kContinuous)) continue;
         procStep = td->GetProcStep(iproc);
         pstep = procStep[i];
         if (pstep < tracks.fPstepV[i]) {
            tracks.fPstepV[i] = pstep;
            tracks.fProcessV[i] = iproc;
         }
      }
   }      
*/
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
   if (!fApplication) {
      Printf("No user application attached - aborting");
      return;
   }   
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
   fWMgr->CreateBaskets();   // geometry should be created by now

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
