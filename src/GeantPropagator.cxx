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
#include "GeantParticleBuffer.h"

GeantPropagator *gPropagator = 0;
   
ClassImp(GeantPropagator)

GeantPropagator *GeantPropagator::fgInstance = 0;

//______________________________________________________________________________
GeantPropagator::GeantPropagator()
                :TObject(),
                 fNthreads(1),
                 fNevents(100),
                 fNtracks(0),
                 fNtransported(0),
                 fNsafeSteps(0),
                 fNsnextSteps(0),
                 fNprocesses(3),
                 fElossInd(0),
                 fNstart(0),
                 fMaxTracks(100000),
                 fMaxThreads(100),
                 fNminThreshold(10),
                 fDebugTrk(-1),
                 fMaxSteps(10000),
                 fNaverage(0.),
                 fEmin(0.1), // 100 MeV
                 fEmax(10),  // 10 Gev
                 fBmag(1.),
                 fUsePhysics(kTRUE),
                 fUseDebug(kFALSE),
                 fTransportOngoing(kFALSE),
                 fSingleTrack(kFALSE),
                 fFillTree(kFALSE),
                 fWMgr(0),
                 fOutput(0),
                 fKineTF1(0),
                 fOutTree(0),
                 fOutFile(0),
                 fTimer(0),
                 fMatrix(0),
                 fVolume(0),
                 fProcesses(0),
                 fTracks(0),
                 fDblArray(0),
                 fProcStep(0),
                 fRndm(0),
                 fPartInd(0),
                 fPartNext(0),
                 fPartTodo(0),
                 fPartCross(0),
                 fFieldPropagator(0),
                 fRotation(0)
{
// Constructor
   fVertex[0] = fVertex[1] = fVertex[2] = 0.;
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
   if (fMatrix) {
     for (i=0; i<fMaxThreads; i++) delete fMatrix[i];
     delete fMatrix;
   }
   delete [] fVolume;
   delete [] fTracks;
   delete [] fDblArray;
   delete [] fProcStep;
   if (fRndm) {
      for (i=0; i<fNthreads; i++) delete fRndm[i];
      delete [] fRndm;
   }
   if (fPartInd) {
      for (i=0; i<fNthreads; i++) delete fPartInd[i];
      delete [] fPartInd;
   }   
   if (fPartNext) {
      for (i=0; i<fNthreads; i++) delete fPartNext[i];
      delete [] fPartNext;
   }   
   if (fPartTodo) {
      for (i=0; i<fNthreads; i++) delete fPartTodo[i];
      delete [] fPartTodo;
   }   
   if (fPartCross) {
      for (i=0; i<fNthreads; i++) delete fPartCross[i];
      delete [] fPartCross;
   }  
   if (fFieldPropagator) {
      for (i=0; i<fNthreads; i++) delete fFieldPropagator[i];
      delete fFieldPropagator;
   }
   if (fRotation) {
      for (i=0; i<fNthreads; i++) delete fRotation[i];
      delete fRotation;
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
   TThread::Lock();
   Int_t iret;
   track->particle = fNtracks;
   fTracks[fNtracks] = track;
//   Int_t tid = TGeoManager::ThreadId();
   fNtracks++;
   fNtransported++;
   if (fNtracks==fMaxTracks) {
      GeantTrack **array = new GeantTrack*[2*fMaxTracks];
      memcpy(array, fTracks, fNtracks*sizeof(GeantTrack*));
      delete [] fTracks;
      fTracks = array;
      fMaxTracks *= 2;
      // Other arrays will need to be also increased...
   }
   iret = fNtracks-1;
   TThread::UnLock();
   return iret;   
}

//______________________________________________________________________________
GeantVolumeBasket *GeantPropagator::ImportTracks(Int_t nevents, Double_t average)
{
// Import tracks from "somewhere". Here we just generate nevents.
   Int_t tid = TGeoManager::ThreadId();
   TGeoNode *node = gGeoManager->FindNode(fVertex[0], fVertex[1], fVertex[2]);
   *fMatrix[tid] = gGeoManager->GetCurrentMatrix();
   fVolume[tid] = node->GetVolume();
   GeantVolumeBasket *basket = (GeantVolumeBasket*)fVolume[tid]->GetField();
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
   for (Int_t event=0; event<nevents; event++) {
      ntracks = fRndm[tid]->Poisson(average);
      
      for (Int_t i=0; i<ntracks; i++) {
//         TGeoBranchArray *a = new TGeoBranchArray();
//         a->InitFromNavigator(gGeoManager->GetCurrentNavigator());
         GeantTrack *track = new GeantTrack(0);
         *track->path = a;
         track->event = event;
         track->particle = fNstart;
         Double_t prob=fRndm[tid]->Uniform(0.,pdgProb[kMaxPart-1]);
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
         Double_t eta = fRndm[tid]->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*fRndm[tid]->Rndm();
         track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
         track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
         track->pz = p*TMath::Cos(theta);
         track->frombdr = kFALSE;
         AddTrack(track);
         basket->AddTrack(fNstart);
         fNstart++;
      }
      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      for (Int_t i=0; i<kMaxPart; i++)
         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
   }      
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

   if (!fMatrix) {
      fMatrix = new TGeoHMatrix*[fMaxThreads];
      for (Int_t i=0; i<fMaxThreads; i++) fMatrix[i] = new TGeoHMatrix();
   } 
   if (!fVolume) {
      fVolume = new TGeoVolume*[fMaxThreads];
      memset(fVolume, 0, fMaxThreads*sizeof(TGeoVolume*));
   }   
   if (!fTracks) {
      fTracks = new GeantTrack*[fMaxTracks];
      memset(fTracks, 0, fMaxTracks*sizeof(GeantTrack*));
   }   
   if (!fDblArray) {
      fDblArray = new Double_t[5*fMaxTracks];
      memset(fDblArray, 0, 5*fMaxTracks*sizeof(Double_t));
   }   
   if (!fProcStep) {
      fProcStep = new Double_t[fNprocesses*fMaxTracks];
      memset(fProcStep, 0, fNprocesses*fMaxTracks*sizeof(Double_t));
   }   
   if (!fRndm) {
      fRndm = new TRandom*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) fRndm[i] = new TRandom();
   }   
   if (!fPartInd) {
      fPartInd = new TArrayI*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) fPartInd[i]   = new TArrayI(fMaxTracks/fNthreads);
   }   
   if (!fPartNext) {
      fPartNext = new TArrayI*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) fPartNext[i]   = new TArrayI(fMaxTracks/fNthreads);
   }   
   if (!fPartTodo) {
      fPartTodo = new TArrayI*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) fPartTodo[i]   = new TArrayI(fMaxTracks/fNthreads);
   }   
   if (!fPartCross) {
      fPartCross = new TArrayI*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) fPartCross[i]   = new TArrayI(fMaxTracks/fNthreads);
   } 
   if (!fRotation) {
      fRotation = new TGeoRotation*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) {
         fRotation[i] = new TGeoRotation();
      }   
   }
   if (!fFieldPropagator) {
      fFieldPropagator = new TGeoHelix*[fNthreads];
      for (Int_t i=0; i<fNthreads; i++) {
         fFieldPropagator[i] = new TGeoHelix(1,1);
         fFieldPropagator[i]->SetField(0,0,fBmag, kFALSE);
      }
   }      
   fWMgr = WorkloadManager::Instance(fNthreads);
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
   // Fill interaction lengths for all processes and all particles
   for (iproc=0; iproc<fNprocesses; iproc++) 
      fProcesses[iproc]->ComputeIntLen(fVolume[tid], ntracks, trackin, &fProcStep[(tid*fNprocesses+iproc)*ntracks],tid);
   // Loop tracks and select process
   for (Int_t i=0; i<ntracks; i++) {
      ipart = trackin[i];
      track = fTracks[ipart];
      track->step = maxlen;
      track->process = -1;
      for (iproc=0; iproc<fNprocesses; iproc++) {
         pstep = fProcStep[(tid*fNprocesses+iproc)*ntracks+ipart];
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
void GeantPropagator::PropagatorGeom(const char *geomfile, Int_t nthreads, Bool_t graphics, Bool_t single, Double_t vertx, Double_t verty, Double_t vertz)
{
// Propagate fNevents in the volume containing the vertex. 
// Simulate 2 physics processes up to exiting the current volume.
   static Bool_t called=kFALSE;
   fNthreads = nthreads;
   fSingleTrack = single;
   Initialize();
   if (called) {
      Printf("Sorry, you can call this only once per session.");
      return;
   }
   called = kTRUE;
   
   // Initialize vertex
   fVertex[0] = vertx;
   fVertex[1] = verty;
   fVertex[2] = vertz;
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

   TCanvas *c1=0;
   TH1F *hnb=0, *hbaskets=0;
   TPad *pad1=0, *pad2=0;
   if (graphics) {
      c1 = new TCanvas("c1","c1",800,900);
      c1->Divide(1,2);
      pad1 = (TPad*)c1->cd(1);
      hnb = new TH1F("hnb","number of baskets per generation",500,0,500);
      hnb->SetFillColor(kRed);
      hnb->Draw();
      pad2 = (TPad*)c1->cd(2);
      hbaskets = new TH1F("hbaskets","baskets population per generation",100,0,100);
      //hbaskets->GetXaxis()->SetTitle("basket");
      hbaskets->SetFillColor(kBlue);
      hbaskets->Draw();
   }
   
   // Create the main volume basket
   GeantVolumeBasket *basket = ImportTracks(fNevents, fNaverage);
//   fBasketArray[fNbaskets++] = basket;
//   fWMgr->AddBasket(basket);

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
//   gGeoManager->SetMultiThread(kTRUE);
   fTimer = new TStopwatch();
   fTimer->Start();
   while (fTransportOngoing) {
      if (fSingleTrack) {
         basket->TransportSingle();
         break;
      }    
      fTransportOngoing = kFALSE;
      // Select baskets to be transported and start transport
      fWMgr->SelectBaskets();
      // Wait for work do get done
      fWMgr->WaitWorkers();
      // Clear transported baskets
      fWMgr->SetFlushed(kFALSE);
      while (!fWMgr->IsFlushed()) {      
         fWMgr->SetFlushed(fWMgr->GetBuffer()->FlushBaskets());
      }   
      fWMgr->GetBuffer()->Reset();
   }   
   fTimer->Stop();
   if (fFillTree) fOutTree->AutoSave();
   delete fOutFile;
   Double_t rtime = fTimer->RealTime();
   Double_t ctime = fTimer->CpuTime();
   fTimer->Print();
   Double_t speedup = ctime/rtime;
   Double_t efficiency = speedup/nthreads;
   fWMgr->Print();
   fWMgr->JoinThreads();
   const char *geomname=geomfile;
   if(strstr(geomfile,"http://root.cern.ch/files/")) geomname=geomfile+strlen("http://root.cern.ch/files/");
   Printf("=== Transported: %lld,  safety steps: %lld,  snext steps: %lld, RT=%gs, CP=%gs", fNtransported, fNsafeSteps, fNsnextSteps,rtime,ctime);
   Printf("   nthreads=%d  speed-up=%f  efficiency=%f", nthreads, speedup, efficiency);
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
