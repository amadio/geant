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
#include "sync_objects.h"

const Int_t     gNevents     = 1;    // number of events transported
const Double_t  gNaverage    = 10000.;  // average number of tracks per event (Poisson)
const Bool_t gUsePhysics  = kTRUE;    // change this to disable physics
Bool_t       gSingleTrack = kFALSE;   // propagate single track versus vectorized

const Bool_t gFillTree    = kFALSE;   // enable I/O

Int_t        gNstart      = 0;        // cumulated initial number of tracks
Int_t        gMaxTracks   = 100000;
Int_t        gNminThreshold = 10;  // threshold for starting transporting a basket
Int_t        gBasketGeneration = 0;
const Int_t  gNprocesses  = 3;
const Int_t  gMaxThreads  = 100;

const Bool_t gUseDebug   = kFALSE;
const Int_t  gDebugTrk   = -1;
Int_t        gNtracks    = 0;
Int_t        gNthreads   = 1;        // number of threads
const Int_t  gMaxSteps   = 10000;      // max number of steps per track

Double_t *gDblArray      = new Double_t[5*gMaxTracks];
Double_t *gProcStep      = new Double_t[gNprocesses*gMaxTracks];
TRandom  **gRndm         = new TRandom*[gMaxThreads];
TArrayI  **gPartInd      = new TArrayI*[gMaxThreads];
TArrayI  **gPartNext     = new TArrayI*[gMaxThreads];
TArrayI  **gPartTodo     = new TArrayI*[gMaxThreads];
TArrayI  **gPartCross    = new TArrayI*[gMaxThreads];
TTree    *gOutTree       = 0;
TFile    *gOutFile       = 0;
TGeoRotation *gRotation  = 0;
Long64_t    gNsafeSteps  = 0;    // number of safety propagations
Long64_t    gNsnextSteps = 0;    // number of snext propagations
Long64_t    gNtransported = 0;    // number of transported particles

enum Species_t {kHadron, kLepton};
enum ProcessType_t {kContinuous, kDiscrete};
enum TrackStatus_t {kAlive, kKilled, kBoundary};

const char *gProcessName[gNprocesses]   = {"Scattering", "Eloss", "Interaction"};
ProcessType_t gProcessType[gNprocesses] = {kDiscrete, kContinuous, kDiscrete};

const Double_t Bmag = 1.;   // kGauss
const Double_t gTolerance = TGeoShape::Tolerance();
// Conversion factor from kGauss to cm (?)
const Double_t kB2C=-0.299792458e-3;
const Double_t emin = 0.1; // 100 MeV
const Double_t emax = 10.0; // 10  GeV
Double_t gVertex[3] = {0., 0., 0.};

TF1        *gKineTF1;
//TPolyMarker3D **pm = new TPolyMarker3D*[gMaxTracks];
//BlockStart  *bufferStart;
//BlockStart  *bufferStop;
concurrent_queue<int> *feeder_queue;
concurrent_queue<int> *answer_queue;

// Array of pointers to ComputeIntLen and Poststep functions
typedef void (*ComputeIntLenFunc)(Int_t ntracks, Int_t *trackin, Double_t *lengths, Int_t tid);
typedef void (*PostStepFunc)(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t* trackout, Int_t tid);
class GeantVolumeBasket;

ComputeIntLenFunc ComputeIntLen[gNprocesses] = {NULL};
PostStepFunc      PostStepAction[gNprocesses] = {NULL};

// Current volume
TGeoVolume *gVolume = 0;
// Number of volume baskets
Int_t gNbaskets = 0;
// Transport sensor
Bool_t gTransportOngoing = kFALSE;
// Current global matrix for the volume
TGeoMatrix *gMatrix = 0;
// Field propagator
TGeoHelix **gFieldPropagator;
TStopwatch gTimer;
//TBits gUsed(gMaxTracks);
//TBits gPushed(gMaxTracks);
//TBits **gUsed;
//TBits **gPushed;

//______________________________________________________________________________
struct GeantTrack {
   Int_t    event;     // event number
   Int_t    particle;  // index of corresponding particle
   Int_t    pdg;       // particle pdg code
   Species_t species;  // particle species
   TrackStatus_t status; // track status
   Int_t    charge;    // particle charge
   Double_t mass;      // particle mass
   Int_t    process;   // current process
   Double_t xpos;      // position
   Double_t ypos;
   Double_t zpos;
   Double_t px;        // momentum
   Double_t py;
   Double_t pz;
   Double_t e;         // energy
   Double_t pstep;     // selected physical step
   Double_t step;      // current step
   Double_t snext;     // straight distance to next boundary
   Double_t safety;    // safe distance to any boundary
   Bool_t   frombdr;   // true if starting from boundary
   Int_t    izero;     // number of small steps used to catch errors
   Int_t    nsteps;    // number of steps made
   TGeoBranchArray *path; // path for this particle in the geometry
//   Int_t    pathindex; // Index of the branch array object in the owner basket
   
   GeantTrack() : event(-1),particle(-1),pdg(0),species(kHadron),status(kAlive),charge(0),mass(0),process(-1),xpos(0),ypos(0),zpos(0),px(0),py(0),pz(0),e(0), pstep(1.E20), step(0), snext(0), safety(0), frombdr(false), izero(0), nsteps(0), path(0) {}
   Double_t           Curvature() {return TMath::Abs(kB2C*Bmag/Pt());}
   void               Direction(Double_t dir[3]);
   Bool_t             IsAlive() const {return (status != kKilled);}
   Bool_t             IsOnBoundary() const {return (status == kBoundary);}
   void               Kill()        {status = kKilled;}
   void               Print(Int_t trackindex=0) const;
   Bool_t             PropagateInFieldSingle(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateInField(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateStraight(Double_t step, Int_t itrack);
   Double_t           Pt()    const {return TMath::Sqrt(px*px+py*py);}
   Double_t           P()     const {return TMath::Sqrt(px*px+py*py+pz*pz);}
   Double_t           Gamma() const {return mass?e/mass:TMath::Limits<double>::Max();}
   Double_t           Beta()  const {return P()/e;}
};

// Array of input tracks
GeantTrack **gTracks = new GeantTrack*[gMaxTracks];

//______________________________________________________________________________
Int_t AddTrack(GeantTrack *track)
{
// Add a new track in the system.
   TThread::Lock();
   Int_t iret;
   track->particle = gNtracks;
   gTracks[gNtracks] = track;
//   Int_t tid = TGeoManager::ThreadId();
//   gUsed[tid]->SetBitNumber(gNtracks, kFALSE);
   gNtracks++;
   gNtransported++;
   if (gNtracks==gMaxTracks) {
      GeantTrack **array = new GeantTrack*[2*gMaxTracks];
      memcpy(array, gTracks, gNtracks*sizeof(GeantTrack*));
      delete [] gTracks;
      gTracks = array;
      gMaxTracks *= 2;
      // Other arrays will need to be also increased...
   }
   iret = gNtracks-1;
   TThread::UnLock();
   return iret;   
}

//______________________________________________________________________________
Double_t BetheBloch(GeantTrack* track, Double_t tz, Double_t ta, Double_t rho) 
{
  if (tz<1. || ta<1.) return 0.;
  const Double_t konst = 0.1535; // MeV cm2/g
  const Double_t emass = 1000*TDatabasePDG::Instance()->GetParticle(kElectron)->Mass();
  const Double_t beta = track->Beta();
  const Double_t gamma = track->Gamma();
  const Double_t bg = beta*gamma;
  const Double_t wmax = 2*emass*bg*bg;
  Double_t ioniz;
  if(tz<13) ioniz = 12 + 7/tz;
  else ioniz = 9.76 + 58.8*TMath::Power(tz,-1.19);

  Double_t bethe = (konst * tz * rho * track->charge * track->charge)/(ta * beta * beta);
//  Printf("ioniz %f",ioniz);
  bethe *= TMath::Log(2*emass*bg*bg*wmax*1e12/(ioniz*ioniz))-2*beta*beta;
//  Printf("bethe %f",bethe);
  return 1.e-3*bethe;
}

//______________________________________________________________________________
Double_t bbf1(Double_t *x, Double_t *par)
{
  Double_t pimass = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  GeantTrack t;
  Double_t bg = TMath::Power(10.,*x);
  Double_t gamma = TMath::Sqrt(bg*bg+1);

  t.px = bg*pimass;
  t.py = 0;
  t.pz = 0;
  t.e = gamma*pimass;
  t.charge = 1;
  t.mass = pimass;
  return 1000*BetheBloch(&t,par[0],par[1],par[2]);
}

//______________________________________________________________________________
void PlotBB(Double_t z, Double_t a, Double_t rho, Double_t bgmin=1e-2, Double_t bgmax=1e6)
{
  TF1 *f=new TF1("bb",bbf1,TMath::Log10(bgmin),TMath::Log10(bgmax),3);
  TH1F *h=new TH1F("hh","Bethe Bloch",100,TMath::Log10(bgmin),TMath::Log10(bgmax));
  h->SetMinimum(1.);
  h->SetMaximum(500.);
  f->SetParameter(0,z);
  f->SetParameter(1,a);
  f->SetParameter(2,rho);
  h->Draw();
  f->Draw("same");
}
//______________________________________________________________________________
void ResetStep(Int_t ntracks, Int_t *array)
{
// Reset current step for a list of tracks.
   for (Int_t i=0; i<ntracks; i++) {
      GeantTrack *track = gTracks[array[i]];
      track->step = 0.;
   }
}

//______________________________________________________________________________
void StepManager(Int_t iproc, Int_t npart, Int_t */*particles*/, Int_t nout, Int_t */*partnext*/)
{
// User stepping routine. <partnext> array can
// be null.
//   for (Int_t ipart=0; ipart<npart; ipart++) gTracks[particles[ipart]]->nsteps++;
   if (gUseDebug) {
      Printf("StepManager: process %s, npart=%d, nout=%d", gProcessName[iproc], npart, nout);
   }   
}

//______________________________________________________________________________
class GeantVolumeBasket : public TObject {
protected:
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNtracks;               // Number of tracks
   Int_t             fFirstFree;             // First un-processed track
   Int_t             fMaxTracks;             // Max number of tracks
   Int_t             fNpaths;                // Max number of physical node paths
   Int_t             fNpathEntries;          // Number of physical paths stored in the array
   Int_t            *fIndex;                 //[fNtracks] Track indices in the global stack
   TGeoBranchArray **fPaths;                 //! Array of physical node paths

   Int_t             InsertPath(TGeoBranchArray* a);
public:
   GeantVolumeBasket(TGeoVolume *vol) : TObject(), fVolume(vol), fNtracks(0), fFirstFree(0), fMaxTracks(50), fNpaths(10), fNpathEntries(0), fIndex(0), fPaths(0) {} 
   virtual ~GeantVolumeBasket();
   
   void              AddTrack(Int_t itrack, TGeoBranchArray* a);
   virtual void      Clear(Option_t *option="");
   void              ComputeTransportLengthSingle(Int_t *trackin);
   void              ComputeTransportLength(Int_t ntracks, Int_t *trackin);
   TGeoBranchArray  *GetBranchArray(GeantTrack *track) const {return track->path;}
   TGeoBranchArray  *GetBranchArray(Int_t itrack) const {return gTracks[fIndex[itrack]]->path;}
   Int_t             GetIndex(Int_t itrack) const {return fIndex[itrack];}
   const Int_t      *GetIndArray() const          {return fIndex;}
   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNtracks() const           {return fNtracks;}
   GeantTrack       *GetTrack(Int_t itrack) const {return gTracks[fIndex[itrack]];}
   TGeoVolume       *GetVolume() const            {return fVolume;}
   void              GetWorkload(Int_t &indmin, Int_t &indmax);
   virtual void      Print(Option_t *option="") const;
   Bool_t            PropagateTrack(Int_t *trackin);
   void              PropagateTracks(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t *trackout, Int_t &ntodo, Int_t *tracktodo, Int_t &ncross, Int_t *trackcross);
   void              TransportSingle();
//   void              TransportTracks();
   
   ClassDef(GeantVolumeBasket,1)  // A path in geometry represented by the array of indices
};   

// Current volume basket
GeantVolumeBasket *gCurrentBasket = 0;
GeantVolumeBasket ** gBasketArray = 0;

// #######   Output class used for I/O #########################################
//______________________________________________________________________________
class GeantOutput : public TObject {
public:
   Double_t        fCpuTime;                 // Cpu time
   Int_t           fVolId;                   // Volume transporting this generation
   Int_t           fBasketGeneration;        // Burent generation of baskets to be flushed
   Int_t           fGeneration;              // Current generation for one basket
   Int_t           fNtracks;                 // Number of tracks in current generation
   Int_t          *fEvent;                   //[fNtracks]
   Int_t          *fInd;                     //[fNtracks] Track indices
   Int_t          *fProc;                    //[fNtracks] Selected processes for each track
   Double_t       *fX;                       //[fNtracks] X positions
   Double_t       *fY;                       //[fNtracks] Y positions
   Double_t       *fZ;                       //[fNtracks] Z positions
   Double_t       *fPx;                      //[fNtracks] Px
   Double_t       *fPy;                      //[fNtracks] Py
   Double_t       *fPz;                      //[fNtracks] Pz
   Double_t       *fE;                       //[fNtracks] E
   Double_t       *fPstep;                   //[fNtracks] Physics step selected
   Double_t       *fStep;                    //[fNtracks] Current step
   Double_t       *fSnext;                   //[fNtracks] Snext distance
   Double_t       *fSafety;                  //[fNtracks] Snext distance

public:
   GeantOutput() : TObject(),fCpuTime(0),fVolId(-1),fBasketGeneration(0),fGeneration(0),fNtracks(0),fEvent(0),fInd(0),fProc(0),fX(0),fY(0),fZ(0),fPx(0),fPy(0),fPz(0),fE(0),fPstep(0),fStep(0),fSnext(0),fSafety(0) {}
   virtual ~GeantOutput();
   void            Init(Int_t size);
   void            Reset();
   void            SetStamp(Int_t volId, Int_t basket_gen, Int_t generation, Int_t ntracks, Double_t cputime=0.) {fVolId=volId; fBasketGeneration=basket_gen; fGeneration=generation;fNtracks=ntracks;fCpuTime=cputime;}
   void            SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t e, Double_t pstep, Double_t step, Double_t snext, Double_t safety);
   void            SetTrack(Int_t ntrack, GeantTrack *track);
ClassDef(GeantOutput,1)       // The transport output per generation
};

GeantOutput *gOutdata = 0;

ClassImp(GeantOutput)

//______________________________________________________________________________
GeantOutput::~GeantOutput()
{
// Destructor
   Reset();
}

//______________________________________________________________________________
void GeantOutput::Init(Int_t size)
{
// Initialize arrays to a given size.
   Reset();
   fEvent = new Int_t[size];
   fInd = new Int_t[size];
   fProc = new Int_t[size];
   fX = new Double_t[size];
   fY = new Double_t[size];
   fZ = new Double_t[size];
   fPx = new Double_t[size];
   fPy = new Double_t[size];
   fPz = new Double_t[size];
   fE  = new Double_t[size];
   fPstep = new Double_t[size];
   fStep = new Double_t[size];
   fSnext = new Double_t[size];
   fSafety = new Double_t[size];
}   
   
//______________________________________________________________________________
void GeantOutput::Reset()
{
// Reset arrays
   delete [] fEvent; fEvent = 0;
   delete [] fInd; fInd = 0;
   delete [] fProc;
   delete [] fX; delete [] fY; delete [] fZ;
   delete [] fPx; delete [] fPy; delete [] fPz; delete [] fE;
   delete [] fPstep; delete [] fStep; delete [] fSnext; delete [] fSafety;
}   

//______________________________________________________________________________
void GeantOutput::SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t e, Double_t pstep, Double_t step, Double_t snext, Double_t safety)
{
// Set parameters for ntrack
   fInd[ntrack] = itrack;
   fEvent[ntrack] = event;
   fProc[ntrack] = proc;
   fX[ntrack] = x;
   fY[ntrack] = y;
   fZ[ntrack] = z;
   fPx[ntrack] = px;
   fPy[ntrack] = py;
   fPz[ntrack] = pz;
   fE[ntrack] = e;
   fPstep[ntrack] = pstep;
   fStep[ntrack] = step;
   fSnext[ntrack] = snext;
   fSafety[ntrack] = safety;
}   

//______________________________________________________________________________
void GeantOutput::SetTrack(Int_t ntrack, GeantTrack *track)
{
// Set parameters for ntrack based on a GeantTrack
   SetTrack(ntrack, track->particle, track->event, track->process, track->xpos, track->ypos, track->zpos, track->px, track->py, track->pz, track->e, track->pstep, track->step, track->snext, track->safety);
}   

//#### Track methods ########################################

//______________________________________________________________________________
void GeantTrack::Direction(Double_t dir[3]) {
   dir[0] = px; dir[1] = py; dir[2] = pz;
   TMath::Normalize(dir);
}

//______________________________________________________________________________
void GeantTrack::Print(Int_t) const {
   TString spath;
//   if (path) path->GetPath(spath);
   Printf("=== Track %d (%s): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Mom:(%f,%f,%f) P:%g E:%g snext=%g safety=%g nsteps=%d",
           particle,spath.Data(), process,pstep,charge,xpos,ypos,zpos,px,py,pz,TMath::Sqrt(px*px+py*py+pz*pz),e,snext,safety,nsteps);
}

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateStraight(Double_t crtstep, Int_t itr)
{
// Propagate along a straight line for neutral particles, for B=0 or for last tiny step.
// The method adds the particle to the next volume basket. 
// Returns the basket pointer, null if exiting geometry.
// Find next volume
   static Int_t istep;
   istep++;
   Double_t dir[3];
   frombdr = kTRUE;
   Direction(dir);
   pstep -= crtstep;
   safety = 0;
   // Change path to reflect the physical volume for the current track;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
//   Int_t tid = nav->GetThreadId();
   TGeoBranchArray *a = path;
   a->UpdateNavigator(nav);
   nav->SetOutside(kFALSE);
   nav->SetStep(crtstep);
   xpos += crtstep*dir[0];
   ypos += crtstep*dir[1];
   zpos += crtstep*dir[2];
   nav->SetCurrentPoint(xpos,ypos,zpos);
   nav->SetCurrentDirection(dir);
   TGeoNode *skip = nav->GetCurrentNode();
   TGeoNode *next = 0;
   next = nav->CrossBoundaryAndLocate(kTRUE, skip);
   if (!next) {
      a->UpdateNavigator(nav);
      next = nav->FindNextBoundaryAndStep(TGeoShape::Big(),kFALSE);
   }   
   gGeoManager->SetVerboseLevel(0);
   xpos += 10.*gTolerance*dir[0];
   ypos += 10.*gTolerance*dir[1];
   zpos += 10.*gTolerance*dir[2];
   if (nav->IsOutside()) return 0;
   // Create a new branch array
   TGeoBranchArray *b = new TGeoBranchArray();
   b->InitFromNavigator(nav);
   TGeoVolume *vol = nav->GetCurrentVolume();
   if (vol->IsAssembly()) Printf("### ERROR ### Entered assembly %s", vol->GetName());
   GeantVolumeBasket *basket = 0;
   if (!vol->GetField()) {
      basket = new GeantVolumeBasket(vol);
      vol->SetField(basket);
      gBasketArray[gNbaskets++] = basket;
   } 
   basket = (GeantVolumeBasket*)vol->GetField();
   basket->AddTrack(itr, b);
   // Signal that the transport is still ongoing if the particle entered a new basket
   if (basket!=gCurrentBasket) {
      gTransportOngoing = kTRUE;
//      gPushed[tid]->SetBitNumber(itr,kTRUE);
   }   
   if (gUseDebug && (gDebugTrk==itr || gDebugTrk<0)) {
      Printf("   track %d: entering %s at:(%f, %f, %f)", itr, nav->GetCurrentNode()->GetName(), xpos,ypos,zpos);
      if (gTracks[itr]->path) gTracks[itr]->path->Print();
   }
//   basket->Print();
   return basket;  
}   

//______________________________________________________________________________
GeantVolumeBasket *GeantTrack::PropagateInField(Double_t crtstep, Bool_t checkcross, Int_t itr)
{
// Propagate with step using the helix propagator. Returns a basket if a 
// boundary was crossed. In such case, the track position and step will reflect
// the boundary crossing point.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Int_t tid = nav->GetThreadId();
   nav->ResetState();
   if (checkcross) {
      path->UpdateNavigator(nav);
      nav->SetLastSafetyForPoint(safety, &xpos);
   }   
   // Reset relevant variables
   frombdr = kFALSE;
   pstep -= crtstep;
   safety -= crtstep;
   if (safety<0.) safety = 0.;
   step += crtstep;
   // Set curvature, charge
   Double_t c = Curvature();
   gFieldPropagator[tid]->SetXYcurvature(c);
   gFieldPropagator[tid]->SetCharge(charge);
   gFieldPropagator[tid]->SetHelixStep(TMath::Abs(TMath::TwoPi()*pz/(c*Pt())));
   gFieldPropagator[tid]->InitPoint(xpos,ypos,zpos);
   Double_t dir[3];
   Direction(dir);
   gFieldPropagator[tid]->InitDirection(dir);
   gFieldPropagator[tid]->UpdateHelix();
   gFieldPropagator[tid]->Step(crtstep);
   const Double_t *point = gFieldPropagator[tid]->GetCurrentPoint();
   const Double_t *newdir = gFieldPropagator[tid]->GetCurrentDirection();
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   Double_t ptot = P();
   px = ptot*newdir[0];
   py = ptot*newdir[1];
   pz = ptot*newdir[2];
   if (!checkcross) return 0;
   if (nav->IsSameLocation(xpos,ypos,zpos,kTRUE)) return 0;
   // Boundary crossed
   TGeoNode *checked = nav->GetCurrentNode();
   TGeoVolume *vol = checked->GetVolume();
   Double_t ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;
   Bool_t outside = nav->IsOutside();
   // Swap track direction and compute distance back to boundary
   dir[0] = -newdir[0]; dir[1] = -newdir[1]; dir[2] = -newdir[2];
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   TGeoNode *node1 = 0;
   TGeoNode *node2 = 0;
   if (level < path->GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         node1 = nav->GetMother(level-lev);
         node2 = path->GetNode(lev);
         if (node1 == node2) {
            if (lev==level) entering = kFALSE;
         } else {
         // different nodes at some level -> entering current node
            break;
         }   
      }
   }   
   if (!entering) checked = path->GetNode(level+1);
   nav->MasterToLocal(&xpos, local);
   nav->MasterToLocalVect(dir, ldir);
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }   
   if (gUseDebug && (gDebugTrk<0 || itr==gDebugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }   
   if (delta>crtstep) {
      if (gUseDebug && (gDebugTrk<0 || itr==gDebugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }   
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;
   // Propagate back to boundary and store position/direction
   xpos += delta*dir[0];
   ypos += delta*dir[1];
   zpos += delta*dir[2]; 
   px = -ptot*dir[0];
   py = -ptot*dir[1];
   pz = -ptot*dir[2];
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   // Mark track as "on boundary" and update step/pstep
   frombdr = kTRUE;
   safety = 0.;
   pstep += delta;
   step -= delta;
   if (outside) {
      nav->SetOutside(kTRUE);
      return 0;
   }   
   // Create a new branch array
   TGeoBranchArray *b = new TGeoBranchArray();
   b->InitFromNavigator(nav);
   if (vol->IsAssembly()) Printf("### ERROR ### Entered assembly %s", vol->GetName());
   GeantVolumeBasket *basket = 0;
   if (!vol->GetField()) {
      basket = new GeantVolumeBasket(vol);
      vol->SetField(basket);
      gBasketArray[gNbaskets++] = basket;
   } 
   basket = (GeantVolumeBasket*)vol->GetField();
   basket->AddTrack(itr, b);
   // Signal that the transport is still ongoing if the particle entered a new basket
   if (basket!=gCurrentBasket) gTransportOngoing = kTRUE;
//   if (gUseDebug && (gDebugTrk==itr || gDebugTrk<0)) {
//      Printf("   track %d: entering %s at:(%f, %f, %f)   ", itr, nav->GetCurrentNode()->GetName(), xpos,ypos,zpos);
//      basket->GetBranchArray(gTracks[itr])->Print();
//   }
   return basket;

}   
//______________________________________________________________________________
Bool_t GeantTrack::PropagateInFieldSingle(Double_t crtstep, Bool_t checkcross, Int_t itr)
{
// Propagate with step using the helix propagator. Navigation must point to
// current particle location. Returns true if crossed next boundary.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Int_t tid = nav->GetThreadId();
   nav->ResetState();
   TGeoBranchArray a;
   a.InitFromNavigator(nav);
   if (checkcross) {
      nav->SetLastSafetyForPoint(safety, &xpos);
   }   
   // Reset relevant variables
   frombdr = kFALSE;
   pstep -= crtstep;
   safety -= crtstep;
   if (safety<0.) safety = 0.;
   step += crtstep;
   // Set curvature, charge
   Double_t c = Curvature();
   gFieldPropagator[tid]->SetXYcurvature(c);
   gFieldPropagator[tid]->SetCharge(charge);
   gFieldPropagator[tid]->SetHelixStep(TMath::Abs(TMath::TwoPi()*pz/(c*Pt())));
   gFieldPropagator[tid]->InitPoint(xpos,ypos,zpos);
   Double_t dir[3];
   Direction(dir);
   gFieldPropagator[tid]->InitDirection(dir);
   gFieldPropagator[tid]->UpdateHelix();
   gFieldPropagator[tid]->Step(crtstep);
   const Double_t *point = gFieldPropagator[tid]->GetCurrentPoint();
   const Double_t *newdir = gFieldPropagator[tid]->GetCurrentDirection();
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   Double_t ptot = P();
   px = ptot*newdir[0];
   py = ptot*newdir[1];
   pz = ptot*newdir[2];
   if (!checkcross) return kFALSE;
   if (nav->IsSameLocation(xpos,ypos,zpos,kTRUE)) return kFALSE;
   // Boundary crossed
   TGeoNode *checked = nav->GetCurrentNode();
   TGeoVolume *vol = checked->GetVolume();
   Double_t ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;
   Bool_t outside = nav->IsOutside();
   // Swap track direction and compute distance back to boundary
   dir[0] = -newdir[0]; dir[1] = -newdir[1]; dir[2] = -newdir[2];
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   TGeoNode *node1 = 0;
   TGeoNode *node2 = 0;
   if (level < a.GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         node1 = nav->GetMother(level-lev);
         node2 = a.GetNode(lev);
         if (node1 == node2) {
            if (lev==level) entering = kFALSE;
         } else {
         // different nodes at some level -> entering current node
            break;
         }   
      }
   }
   if (!entering) checked = a.GetNode(level+1);
   nav->MasterToLocal(&xpos, local);
   nav->MasterToLocalVect(dir, ldir);
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }   
   if (gUseDebug && (gDebugTrk<0 || itr==gDebugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }   
   if (delta>crtstep) {
      if (gUseDebug && (gDebugTrk<0 || itr==gDebugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }   
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;
   // Propagate back to boundary and store position/direction
   xpos += delta*dir[0];
   ypos += delta*dir[1];
   zpos += delta*dir[2]; 
   px = -ptot*dir[0];
   py = -ptot*dir[1];
   pz = -ptot*dir[2];
   xpos = point[0]; ypos = point[1]; zpos = point[2];
   // Mark track as "on boundary" and update step/pstep
   frombdr = kTRUE;
   safety = 0.;
   pstep += delta;
   step -= delta;
   if (outside) {
      nav->SetOutside(kTRUE);
   }
   // Boundary crossed, navigation points to new location
   return kTRUE;
}   


//######## Initialization of tracks and geometry

//______________________________________________________________________________
GeantVolumeBasket *ImportTracks(Int_t nevents, Double_t average)
{
// Import tracks from "somewhere". Here we just generate nevents.
   TGeoNode *node = gGeoManager->FindNode(gVertex[0], gVertex[1], gVertex[2]);
   gMatrix = gGeoManager->GetCurrentMatrix();
   gVolume = node->GetVolume();
   GeantVolumeBasket *basket = new GeantVolumeBasket(gVolume);
   gVolume->SetField(basket);
   
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
      ntracks = gRandom->Poisson(average);
      
      for (Int_t i=0; i<ntracks; i++) {
         TGeoBranchArray *a = new TGeoBranchArray();
         a->InitFromNavigator(gGeoManager->GetCurrentNavigator());
         GeantTrack *track = new GeantTrack();
         track->event = event;
         track->particle = gNstart;
         Double_t prob=gRandom->Uniform(0.,pdgProb[kMaxPart-1]);
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
         track->xpos = gVertex[0];
         track->ypos = gVertex[1];
         track->zpos = gVertex[2];
         track->e = gKineTF1->GetRandom()+track->mass;
         Double_t p = TMath::Sqrt((track->e-track->mass)*(track->e+track->mass));
         Double_t eta = gRandom->Uniform(etamin,etamax);  //multiplicity is flat in rapidity
         Double_t theta = 2*TMath::ATan(TMath::Exp(-eta));
         //Double_t theta = TMath::ACos((1.-2.*gRandom->Rndm()));
         Double_t phi = TMath::TwoPi()*gRandom->Rndm();
         track->px = p*TMath::Sin(theta)*TMath::Cos(phi);
         track->py = p*TMath::Sin(theta)*TMath::Sin(phi);
         track->pz = p*TMath::Cos(theta);
         track->frombdr = kFALSE;
         AddTrack(track);
         basket->AddTrack(gNstart, a);
         gNstart++;
      }
      Printf("Event #%d: Generated species for %6d particles:", event, ntracks);
      for (Int_t i=0; i<kMaxPart; i++)
         Printf("%15s : %6d particles", TDatabasePDG::Instance()->GetParticle(pdgGen[i])->GetName(), pdgCount[i]);
   }      
   return basket;
}

//______________________________________________________________________________
void PrintParticles(Int_t *trackin, Int_t ntracks, Int_t tid)
{
// Print the detailed particles list.
   Printf("================ THREAD %d: particles list", tid);
   for (Int_t i=0; i<ntracks; i++) {
      gTracks[trackin[i]]->Print();
   }
}
      
//______________________________________________________________________________
Bool_t LoadGeometry(const char *filename="geometry.root")
{
// Load the detector geometry from file.
   if (gGeoManager) return kTRUE;
   TGeoManager *geom = TGeoManager::Import(filename);
   if (geom) {
      return kTRUE;
   }   
   ::Error("LoadGeometry","Cannot load geometry from file %s", filename);
   return kFALSE;
}
   
//####### Generate interaction lengths for a bunch of tracks
//______________________________________________________________________________
void ComputeIntLenScattering(Int_t ntracks, Int_t *trackin, Double_t *lengths, Int_t tid)
{
// Generates an interaction length for the scattering process. Nothing physical,
// just generate something comparable with the size of the current volume.
//
// trackin and lengths arrays should be be correctly set by the caller
   const Double_t kC1 = 500.;
   const Double_t xlen = TMath::Limits<double>::Max();
   Int_t itrack;
   Double_t density = 1.e-5;
   TGeoMaterial *mat = gVolume->GetMaterial();
   if (mat) density = mat->GetDensity();
   density = TMath::Max(density, 1.E-3);
   // Make sure we write in the thread space for the current basket
   Double_t *rndArray = &gDblArray[2*tid*ntracks];
   Int_t irnd = 0;
   gRndm[tid]->RndmArray(ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      if (gTracks[itrack]->IsAlive())
        lengths[itrack] = kC1*gTracks[itrack]->e*rndArray[irnd++]/density;
      else
        lengths[itrack] =  0.5*xlen; 
   }
}

//______________________________________________________________________________
void ComputeIntLenEloss(Int_t ntracks, Int_t *trackin, Double_t *lengths, Int_t)
{
// Energy loss process. Continuous process. Compute step limit for losing
// maximum dw per step.
   const Double_t dw = 1.E-3;  // 1 MEV
   TGeoMaterial *mat = gVolume->GetMaterial();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;
   Int_t itrack;
   GeantTrack *track;
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];  
      track = gTracks[itrack];
      if(track->charge && !invalid_material && track->IsAlive()) {
         Double_t dedx = BetheBloch(track,matz,mata,matr);
         Double_t stepmax = (dedx>1.E-32)?dw/dedx:0.5*TMath::Limits<double>::Max();
         lengths[itrack] = stepmax;
      } else {
         lengths[itrack]=0.5*TMath::Limits<double>::Max();
      }      
   }
}

//______________________________________________________________________________
void ComputeIntLenInteraction(Int_t ntracks, Int_t *trackin, Double_t *lengths, Int_t)
{
// 
   Double_t fact = 1.;
   const Double_t nabarn = fact*TMath::Na()*1e-24;
   Int_t itrack;
   GeantTrack *track;
   Double_t xlen = TMath::Limits<double>::Max();
   TGeoMaterial *mat = gVolume->GetMaterial();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;
   if (!invalid_material) {
      Double_t density = TMath::Max(matr,1e-5);
      Double_t sigma = 28.5*TMath::Power(mata,0.75);
      xlen = mat->GetA()/(sigma*density*nabarn);
   } else {
      for (itrack=0; itrack<ntracks; itrack++) lengths[itrack] = 0.5*TMath::Limits<double>::Max();
      return;
   }   
 
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gTracks[itrack];
      if(track->species == kHadron && track->IsAlive()) {
         Double_t ek = track->e - track->mass;
         lengths[itrack] = xlen*(0.007+0.1*TMath::Log(ek)/ek+0.2/(ek*ek));
      } else {
         lengths[itrack] = 0.5*TMath::Limits<double>::Max();
      }
   }
}

//####### Post step actions for different processes, acting for a bunch of tracks
//______________________________________________________________________________
void PostStepScattering(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t* trackout, Int_t tid)
{
// Do post-step actions on particle after scattering process. Surviving tracks
// copied in trackout
   // Compute the max theta angle opening after scattering.
   const Double_t ctmax = TMath::Cos(1.*TMath::DegToRad()); // 1 degree
   Double_t theta, phi, scale,thetav,phiv; 
   Double_t dir[3];
   Double_t dirnew[3];
   GeantTrack *track = 0;
   Int_t itrack;
   Double_t p;
   Double_t *rndArray = &gDblArray[2*tid*ntracks];
   Int_t irnd = 0;
   gRndm[tid]->RndmArray(2*ntracks, rndArray);
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gTracks[itrack];
      theta = TMath::ACos((1.-rndArray[irnd++]*(1.-ctmax)));
      // Re-scale from emin to emax
      scale = (track->e-emin)/(emax-emin);
      theta *= 1-scale;  // hi-energy don't scatter much
      phi = TMath::TwoPi()*rndArray[irnd++];
      // System along the (px,py,pz)
      p = track->P();
      thetav = TMath::ACos(track->pz/p)*TMath::RadToDeg();
      phiv = TMath::ATan2(track->py,track->px)*TMath::RadToDeg();
      gRotation->SetAngles(phiv-90,-thetav,0);
      dir[0] = TMath::Sin(theta)*TMath::Cos(phi);
      dir[1] = TMath::Sin(theta)*TMath::Sin(phi);
      dir[2] = TMath::Cos(theta);
      gRotation->LocalToMaster(dir, dirnew);
      track->px = p*dirnew[0];
      track->py = p*dirnew[1];
      track->pz = p*dirnew[2];
      // All tracks survive
      if (trackout) trackout[nout] = itrack;
      nout++;
   }   
   StepManager(0, ntracks, trackin, nout, trackout);
}   

//______________________________________________________________________________
void PostStepEloss(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t* trackout, Int_t)
{
// Do post-step actions after energy loss process. 
   Double_t eloss, dedx;
   GeantTrack *track;
   Int_t itrack;
   TGeoMaterial *mat = gVolume->GetMaterial();
   Double_t mata = mat->GetA();
   Double_t matz = mat->GetZ();
   Double_t matr = mat->GetDensity();
   Bool_t invalid_material = kFALSE;
   if (matz<1 || mata<1 || matr<1.E-8) invalid_material = kTRUE;

   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];   
      track = gTracks[itrack];
      if (!track->IsAlive()) continue;
      if (track->e-track->mass < emin) {
         track->Kill();
         continue;
      }   
      if (track->step==0 || invalid_material) {
         if (trackout) trackout[nout] = itrack;
         nout++;
         continue;
      }   
      dedx = BetheBloch(track,matz,mata,matr);
      eloss = track->step*dedx;
      if (track->e-track->mass-eloss < emin) eloss = track->e-track->mass;
      Double_t gammaold = track->Gamma();
      Double_t bgold = TMath::Sqrt((gammaold-1)*(gammaold+1));
      track->e -= eloss;
      if (track->e-track->mass < emin) {
         track->Kill();
         continue;
      }   
      if (trackout) trackout[nout] = itrack;
      nout++;

      Double_t gammanew = track->Gamma();
      Double_t bgnew = TMath::Sqrt((gammanew-1)*(gammanew+1));
      Double_t pnorm = bgnew/bgold;
      track->px *= pnorm;
      track->py *= pnorm;
      track->pz *= pnorm;
   }   
   StepManager(1, ntracks, trackin, nout, trackout);
}   

//______________________________________________________________________________
void PostStepInteraction(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t* trackout, Int_t tid)
{
// Do post-step actions on particle after interaction process. 
//   if (gUseDebug) Printf("PostStepInteraction %d tracks", ntracks);
// We calculate the CMS energy
// We suppose at first that the available energy is the Kin cms energy
// We produce equal number of pos and neg pions

   static TGenPhaseSpace gps;
   GeantTrack *track;
   Int_t itrack;
   Double_t *rndArray = &gDblArray[2*tid*ntracks];
   const Double_t pimass = TDatabasePDG::Instance()->GetParticle(kPiMinus)->Mass();
   const Double_t prodm[18] = {pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass,
			       pimass, pimass, pimass, pimass, pimass, pimass};
   gRndm[tid]->RndmArray(ntracks, rndArray);

   Int_t nprod = 0;
   Int_t ngen  = 0;
   for (Int_t i=0; i<ntracks; i++) {
      itrack = trackin[i];
      track = gTracks[itrack];
      Double_t en = track->e;
      Double_t m1 = track->mass;
      Double_t m2 = gVolume->GetMaterial()->GetA();
      Double_t cmsen = TMath::Sqrt(m1*m1+m2*m2+2*en*m2)-m1-m2;
      // Calculate the number of pions as a poisson distribution leaving half of the cms energy
      // for phase space momentum
      Int_t npi = 0.5*gRndm[tid]->Rndm()*cmsen/pimass+0.5;
      if(npi>1) {
         do { nprod = TMath::Min(gRndm[tid]->Poisson(npi),9); } 
         while(nprod*pimass*2>cmsen || nprod==0);
//         Printf("Inc en = %f, cms en = %f produced pis = %d",en,cmsen,nprod);
         TLorentzVector pcms(track->px, track->py, track->pz, track->e + m2);
         if(!gps.SetDecay(pcms,2*nprod,prodm)) Printf("Forbidden decay!");
         gps.Generate();
         //Double_t pxtot=track->px;
         //Double_t pytot=track->py;
         //Double_t pztot=track->pz;
         for(Int_t j=0; j<2*nprod; ++j) {
            TGeoBranchArray *a = new TGeoBranchArray(*track->path);
            GeantTrack *trackg=new GeantTrack();
            TLorentzVector *lv = gps.GetDecay(j);
            if(j%2) trackg->pdg = kPiMinus;
            else trackg->pdg = kPiPlus;
            trackg->species = kHadron;
            trackg->charge = TDatabasePDG::Instance()->GetParticle(trackg->pdg)->Charge()/3.;
            trackg->mass = pimass;
            trackg->process = 0;
            trackg->xpos = track->xpos;
            trackg->ypos = track->ypos;
            trackg->zpos = track->zpos;
            trackg->px = lv->Px();
            trackg->py = lv->Py();
            trackg->pz = lv->Pz();
            trackg->e = lv->E();
//            Double_t mm2 = trackg->e*trackg->e-trackg->px*trackg->px-trackg->py*trackg->py-trackg->pz*trackg->pz;
            Int_t itracknew = AddTrack(trackg);
            trackout[nout++] = trackg->particle = itracknew;
            ngen++;
            if (gCurrentBasket) gCurrentBasket->AddTrack(itracknew,a);
           //check
           //pxtot -= trackg->px;
           //pytot -= trackg->py;
           //pztot -= trackg->pz;
         }
	//	Printf("pbal = %f %f %f",pxtot, pytot, pztot);
      }
   }   
   StepManager(2, ntracks, trackin, nout, trackout);
   if (ngen) {
      // Generated particles may be below threshold-> Call PostStepEloss
      Int_t nsurv = 0;
      Int_t *trackgen = new Int_t[ngen];
      PostStepEloss(ngen, &trackout[nout-ngen], nsurv, trackgen,tid);
      memcpy(&trackout[nout-ngen], trackgen, nsurv*sizeof(Int_t));
      nout += nsurv-ngen;
      delete [] trackgen;
   }
}   

// ########## Selection of the physics process that generated the smallest intlen
//______________________________________________________________________________
void PhysicsSelect(Int_t ntracks, Int_t *trackin, Int_t tid)
{
// Generate all physics steps for the tracks in trackin.
// Vectorized, except the unavoidable Sort()
   static const Double_t maxlen = TMath::Limits<double>::Max();   
   Double_t pstep;
   Int_t ipart, iproc;
   GeantTrack *track;
   // Fill interaction lengths for all processes and all particles
   for (iproc=0; iproc<gNprocesses; iproc++) 
      ComputeIntLen[iproc](ntracks, trackin, &gProcStep[(tid*gNprocesses+iproc)*ntracks],tid);
   // Loop tracks and select process
   for (Int_t i=0; i<ntracks; i++) {
      ipart = trackin[i];
      track = gTracks[ipart];
      track->step = maxlen;
      track->process = -1;
      for (iproc=0; iproc<gNprocesses; iproc++) {
         pstep = gProcStep[(tid*gNprocesses+iproc)*ntracks+ipart];
         if (pstep < track->step) {
            track->step = pstep;
            track->process = iproc;
         }
      }
      if (gUseDebug && (gDebugTrk==ipart || gDebugTrk<0)) {
         Printf("   (%d) PhysicsSelect: track #%d - process=%d pstep=%g",tid,ipart,track->process,track->step);
         if (track->step>1.E200) {
            Printf("xxx");
         }   
      }
   }      
}

// ####### Utilities
//______________________________________________________________________________
void InitFunctions()
{
// Initialize arrays of pointers to functions per process.
   ComputeIntLen[0] = &ComputeIntLenScattering;
   ComputeIntLen[1] = &ComputeIntLenEloss;
   ComputeIntLen[2] = &ComputeIntLenInteraction;
   PostStepAction[0] = &PostStepScattering;
   PostStepAction[1] = &PostStepEloss;
   PostStepAction[2] = &PostStepInteraction;
   gRotation = new TGeoRotation();
   gFieldPropagator = new TGeoHelix*[gNthreads];
//   gUsed   = new TBits*[gNthreads];
//   gPushed = new TBits*[gNthreads];
//   bufferStart = new BlockStart(gNthreads);
//   bufferStop  = new BlockStart(gNthreads);
   feeder_queue = new concurrent_queue<int>(true);
   answer_queue = new concurrent_queue<int>;
   for (Int_t i=0; i<gNthreads; i++) {
      gPartInd[i]   = new TArrayI(gMaxTracks/gNthreads);
      gPartNext[i]  = new TArrayI(gMaxTracks/gNthreads);
      gPartTodo[i]  = new TArrayI(gMaxTracks/gNthreads);
      gPartCross[i] = new TArrayI(gMaxTracks/gNthreads);
      gFieldPropagator[i] = new TGeoHelix(1,1);
      gFieldPropagator[i]->SetField(0,0,Bmag, kFALSE);
      gRndm[i]      = new TRandom();
//      gUsed[i]   = new TBits(gMaxTracks);
//      gPushed[i] = new TBits(gMaxTracks);
   }   
}

//______________________________________________________________________________
void SelectTracksForProcess(Int_t iproc, Int_t ntotransport, Int_t *particles, Int_t &ntodo, Int_t *parttodo)
{
// Add to output array all particles that will do the process iproc.
   for (Int_t itr=0; itr<ntotransport; itr++)
      if (gTracks[particles[itr]]->process == iproc) parttodo[ntodo++] = particles[itr];
}

//#### Volume basket ########################################
ClassImp(GeantVolumeBasket)

//______________________________________________________________________________
GeantVolumeBasket::~GeantVolumeBasket()
{
// Clean up
   delete [] fIndex;
   for (Int_t i=0; i<fNpathEntries; i++) delete fPaths[i];
   delete [] fPaths;
}   

//______________________________________________________________________________
void GeantVolumeBasket::GetWorkload(Int_t &indmin, Int_t &indmax)
{
// Get a range of indices to work with for a given thread.
   TThread::Lock();
   indmin = indmax = fFirstFree;
   if (!fNtracks || fFirstFree==fNtracks) {
      TThread::UnLock();
      return;
   }   
   Int_t fair_share = fNtracks/gNthreads;
   Int_t remaining = fNtracks%gNthreads;
   indmax = indmin+fair_share;
   if (remaining) indmax++;
   if (indmax > fNtracks) indmax = fNtracks;
   fFirstFree = indmax;
   TThread::UnLock();
}   

//______________________________________________________________________________
Int_t GeantVolumeBasket::InsertPath(TGeoBranchArray* a)
{
// Check if the path is stored and return its index. If not there, insert it in fPaths.
// Deletes input branch array if an equal one is found.
   Int_t i = 0;
//   Int_t tid = TGeoManager::ThreadId();
   if (a->GetNode(a->GetLevel())->GetVolume() != fVolume) {
      Printf("ERROR: Wrong path for basket %s", GetName());
   }
   Int_t ifound = TGeoBranchArray::BinarySearch(fNpathEntries,(const TGeoBranchArray**)fPaths,a);
   if (ifound<=fNpathEntries-1) {
      if (ifound>=0 && *a == *fPaths[ifound]) {
         delete a;
         return ifound;
      }   
      // Insert path at ifound+1. Copy elements from ifound+1 up 1 unit
//      memmove(&fPaths[ifound+2], &fPaths[ifound+1], (fNpathEntries-ifound-1)*sizeof(TGeoBranchArray *));
      for (i=fNpathEntries-1; i>ifound; i--) fPaths[i+1]=fPaths[i];
   }
   fPaths[ifound+1] = a;

   // Increase array of paths if needed.
   fNpathEntries++;
   if (fNpathEntries == fNpaths) {
      TGeoBranchArray **newarr = new TGeoBranchArray*[2*fNpaths];
      memcpy(newarr, fPaths, fNpathEntries*sizeof(TGeoBranchArray *));
      fNpaths *=2;
      delete [] fPaths;
      fPaths = newarr;
   }   
   return ifound+1;
}   

//______________________________________________________________________________
void GeantVolumeBasket::AddTrack(Int_t itrack, TGeoBranchArray* a)
{
// Add a track and its path to the basket.
   TThread::Lock();
   if (!fNtracks) {
      fPaths = new TGeoBranchArray*[fNpaths];
      fIndex = new Int_t[fMaxTracks];
      fPaths[fNtracks] = a;
      fIndex[fNtracks] = itrack;
      gTracks[itrack]->path = a;
      fNtracks++;
      fNpathEntries++;
      TThread::UnLock();
      return;
   }
//   Int_t tid = TGeoManager::ThreadId();
   // Check if the path is already stored.
   Int_t ipath = InsertPath(a);
   gTracks[itrack]->path = fPaths[ipath];
   if (gTracks[itrack]->path->GetNode(gTracks[itrack]->path->GetLevel())->GetVolume() != fVolume) {
      Printf("ERROR: Wrong path for basket %s", GetName());
      gTracks[itrack]->Print();
   }  
   TString spath; 
   fPaths[ipath]->GetPath(spath);
//   Printf("(%d) Track %d added to basket:%s in path: %s", tid, itrack, GetName(), spath.Data());
   // If the track is already in the basket, do not add it again (!)
   if (gCurrentBasket != this) {
      fIndex[fNtracks] = itrack;
      fNtracks++;
      // Increase arrays of tracks and path indices if needed
      if (fNtracks == fMaxTracks) {
         Int_t *newindex = new Int_t[2*fMaxTracks];
         memcpy(newindex, fIndex, fNtracks*sizeof(Int_t));
         delete [] fIndex;
         fIndex = newindex;
         fMaxTracks *= 2;
      }   
   }
   TThread::UnLock();   
}     

//______________________________________________________________________________
void GeantVolumeBasket::Clear(Option_t *)
{
// Clear all particles and paths 
   TThread::Lock();
   for (Int_t i=0; i<fNpathEntries; i++) delete fPaths[i];
   fNtracks = 0;
   fFirstFree = 0;
   fNpathEntries = 0;
   TThread::UnLock();   
}   

//______________________________________________________________________________
void GeantVolumeBasket::ComputeTransportLength(Int_t ntracks, Int_t *trackin)
{
// Computes snext and safety for an array of tracks. This is the transportation 
// process. Tracks are assumed to be inside gVolume.
   static Int_t icalls = 0;
   Double_t pdir[3];
   Int_t itr;
   Bool_t isOnBoundary = kFALSE;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   nav->SetOutside(kFALSE);
   
   for (itr=0; itr<ntracks; itr++) {
      GeantTrack *track = gTracks[trackin[itr]];
      track->Direction(pdir);
      track->path->UpdateNavigator(nav);
      nav->SetCurrentPoint(&track->xpos);
      nav->SetCurrentDirection(pdir);
      isOnBoundary = track->frombdr;
      Double_t pstep = TMath::Min(1.E20, track->pstep);
      nav->FindNextBoundary(pstep,"",isOnBoundary);
//      gGeoManager->SetVerboseLevel(0);
      track->safety = nav->GetSafeDistance();
//      if (!isOnBoundary && track->safety<gTolerance) {
//         nav->FindNextBoundary(track->pstep,"",isOnBoundary);
//      }
//      track->snext = nav->GetStep();
      track->snext = TMath::Max(gTolerance,nav->GetStep());
      if (gUseDebug && (gDebugTrk==trackin[itr] || gDebugTrk<0)) {
//         Printf("    %d   track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", icalls,trackin[itr], nav->GetPath(), track->snext, track->safety, track->pstep);
         Printf("       track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", trackin[itr], nav->GetPath(), track->snext, track->safety, track->pstep);
         track->Print(trackin[itr]);
      }   
   }
   icalls++;  
}            

//______________________________________________________________________________
void GeantVolumeBasket::ComputeTransportLengthSingle(Int_t *trackin)
{
// Computes snext and safety for a given track. This is the transportation 
// process. Track is assumed to be inside gVolume. Lighter version of normal 
// navigation, using only volume data and ignoring MANY option.
   static Int_t icalls = 0;
   Double_t pdir[3];
   Bool_t isOnBoundary = kFALSE;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   nav->SetOutside(kFALSE);
   GeantTrack *track = gTracks[trackin[0]];
   track->Direction(pdir);
   nav->SetCurrentPoint(&track->xpos);
   nav->SetCurrentDirection(pdir);
   isOnBoundary = track->frombdr;
   if (gUseDebug && (gDebugTrk==trackin[0] || gDebugTrk<0)) {
//         gGeoManager->SetVerboseLevel(10);
   }   
   Double_t pstep = TMath::Min(1.E20, track->pstep);
   nav->FindNextBoundary(pstep,"",isOnBoundary);
   gGeoManager->SetVerboseLevel(0);
   track->safety = nav->GetSafeDistance();
//   if (!isOnBoundary && track->safety<gTolerance) {
//      nav->FindNextBoundary(track->pstep,"",isOnBoundary);
//   }
//      track->snext = nav->GetStep();
   track->snext = TMath::Max(gTolerance,nav->GetStep());
   if (gUseDebug && (gDebugTrk==trackin[0] || gDebugTrk<0)) {
      Printf("       track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", trackin[0], nav->GetPath(), track->snext, track->safety, track->pstep);
//      Printf("    %d   track %d: %s  snext=%19.15f safe=%19.15f pstep=%f", icalls,trackin[0], nav->GetPath(), track->snext, track->safety, track->pstep);
      track->Print(trackin[0]);
   }   
   icalls++;  
}            

//______________________________________________________________________________
void GeantVolumeBasket::PropagateTracks(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t *trackout, Int_t &ntodo, Int_t *tracktodo, Int_t &ncross, Int_t *trackcross)
{
// Propagate the ntracks with their selected steps. If a boundary is
// found in the way, the track is stopped. Nout must be initialized from outside.
//     trackin = array of <ntracks> input tracks
//     trackout = array of <nout> tracks propagated to physics processes
//     tracktodo = array of <ntodo> tracks propagated with a safe step or crossing 
//                 inside the same volume. These have to be propagated  again.
//     trackcross = array of <ncross> tracks that crossed the boundary. For these tracks
//                 the continuous processes have to be applied after propagation
   GeantTrack *track;
   Double_t step, snext, safety, c;
   ntodo = 0;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   GeantVolumeBasket *basket = 0;
//   Printf("===== PropagateTracks: ntracks=%d nout=%d", ntracks, nout);
   for (Int_t itr=0; itr<ntracks; itr++) {
      track = gTracks[trackin[itr]];
      // Skip neutral tracks for the time being (!!!)
      if (!track->charge) continue;
      if (!track->IsAlive()) continue;
      track->nsteps++;
      if (track->nsteps > gMaxSteps) {
         track->Kill();
         continue;
      }   
      step = track->pstep;
      snext = track->snext;
      safety = track->safety;
      c = track->Curvature();
      // If proposed step less than safety, just propagate to physics process
      if (step<safety) {
         track->izero = 0;
         track->PropagateInField(step, kFALSE, trackin[itr]);
         gNsafeSteps++; // increment-only, thread safe
         // track transported to physics process
         if (gUseDebug && (gDebugTrk==trackin[itr] || gDebugTrk<0)) Printf("   track %d process: %s", trackin[itr], gProcessName[track->process]);
         trackout[nout++] = trackin[itr]; // <- survives geometry, stopped due to physics
         continue; // -> to next track
      }
      // Check if we can propagate to boundary
      if (0.25*c*snext<1E-6 && snext<1E-3 && snext<step-1E-6) {
         // Propagate with snext and check if we crossed
         //   backup track position and momentum
         if (track->izero>10) snext = 1.E-3;
         basket = track->PropagateInField(snext+10*gTolerance, kTRUE, trackin[itr]);
         if (snext<1.E-6) track->izero++;
         else track->izero = 0;
         gNsnextSteps++;
         if (!basket) {
            // track exiting
            if (nav->IsOutside()) {
               if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
               trackcross[ncross++] = trackin[itr];
               continue;
            }   
            // these tracks do not cross
//            if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
            tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
            continue; // -> next track
         }
         // These tracks are reaching boundaries
         trackcross[ncross++] = trackin[itr];
//         if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d pushed to boundary of %s", trackin[itr], basket->GetName());
//         basket = track->PropagateStraight(snext, trackin[itr]);
         if (basket==this) {
            // The track entered the same basket -> add it to the todo list
//           if (gUseDebug && (gDebugTrk==trackin[itr] || gDebugTrk<0)) Printf("   track %d entered the same basket %s", trackin[itr], GetName());
           tracktodo[ntodo++] = trackin[itr]; 
         }   
         continue; // -> to next track
      }
      // Track has safety<pstep but next boundary not close enough.
      // We propagate in field with the safety value.
      if (safety<gTolerance) {
         // Track getting away from boundary. Work to be done here
         // In principle we need a safety value for the content of the current volume only
         // This does not behave well on corners...
         // ... so we peek a small value and chech if this crosses, than recompute safety
         safety = 1.E-3;
         track->izero++;
         if (track->izero > 10) safety = 0.5*snext;
         basket = track->PropagateInField(safety, kTRUE, trackin[itr]);
      } else {
         if (track->izero > 10) {
            // Propagate with snext
            basket = track->PropagateInField(snext+10*gTolerance, kTRUE, trackin[itr]);
            track->izero = 0; 
            gNsnextSteps++;
            if (!basket) {
               // track exiting geometry
               if (nav->IsOutside()) {
                  if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
                  trackcross[ncross++] = trackin[itr];
                  continue;
               }   
               // these tracks do not cross
//               if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
               tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
               continue; // -> next track
            }
            // These tracks are reaching boundaries
            trackcross[ncross++] = trackin[itr];
            if (basket==this) {
               // The track entered the same basket -> add it to the todo list
              tracktodo[ntodo++] = trackin[itr]; 
            }   
            continue;
         }
         if (safety<1.E-3) track->izero++;
         // Propagate with safety without checking crossing
         basket = track->PropagateInField(safety, kFALSE, trackin[itr]);
      }   
      gNsafeSteps++;
      if (!basket) {
         // check if track exiting
         if (nav->IsOutside()) {
            if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
            trackcross[ncross++] = trackin[itr];
            continue;
         }   
         // these tracks do not cross
//         if (gUseDebug && (gDebugTrk== trackin[itr] || gDebugTrk<0)) Printf("   track %d propagated with safety=%19.15f", trackin[itr], safety);
         tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
         continue; // -> next track
      }
      // These tracks are reaching boundaries
      trackcross[ncross++] = trackin[itr];
      if (basket==this) {
         // The track entered the same basket -> add it to the todo list
        if (gUseDebug && (gDebugTrk==trackin[itr] || gDebugTrk<0)) Printf("   track %d entered the same basket %s", trackin[itr], GetName());
        tracktodo[ntodo++] = trackin[itr]; 
      } 
   }   
   // Recompute snext and safety for todo tracks
   if (ntodo) ComputeTransportLength(ntodo, tracktodo);
}

//______________________________________________________________________________
Bool_t GeantVolumeBasket::PropagateTrack(Int_t *trackin)
{
// Propagate a single track with its selected physics step. If a boundary is
// found in the way return. 
   GeantTrack *track;
   Double_t step = 0;
   Double_t snext, safety, c;
//   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   track = gTracks[trackin[0]];
   // Skip neutral tracks for the moment (!!!)
   Bool_t crossed = kFALSE;
   // If proposed step less than safety, just propagate to physics process
   while (!crossed) {
      step = track->pstep;
      snext = track->snext;
      safety = track->safety;
      c = track->Curvature();
      if (step<safety) {
         track->izero = 0;
         crossed = track->PropagateInFieldSingle(step, kFALSE, trackin[0]);
         gNsafeSteps++;
         return crossed;
      }
      // Check if we can propagate to boundary
      if (0.25*c*snext<1E-6 && snext<1E-3 && snext<step-1E-6) {
         // Propagate with snext and check if we crossed
         if (track->izero>10) snext = 1.E-3;
         crossed = track->PropagateInFieldSingle(snext+10*gTolerance, kTRUE, trackin[0]);
         if (snext<1.E-6) track->izero++;
         else track->izero = 0;
         gNsnextSteps++;
         if (crossed) return crossed;
         ComputeTransportLengthSingle(trackin);
         continue;
      }   

      // Track has safety<pstep but next boundary not close enough.
      // We propagate in field with the safety value.
      if (safety<gTolerance) {
         // Track getting away from boundary. Work to be done here
         // In principle we need a safety value for the content of the current volume only
         // This does not behave well on corners...
         // ... so we peek a small value and chech if this crosses, than recompute safety
         safety = 1.E-3;
         track->izero++;
         if (track->izero > 10) safety = 0.5*snext;
         crossed = track->PropagateInFieldSingle(safety, kTRUE, trackin[0]);
      } else {        
         if (track->izero > 10) {
            // Propagate with snext
            crossed = track->PropagateInFieldSingle(snext+10*gTolerance, kTRUE, trackin[0]);
            track->izero = 0; 
            gNsnextSteps++;
            if (crossed) return crossed;
            ComputeTransportLengthSingle(trackin);
            continue;
         }
         if (safety<1.E-3) track->izero++;
         crossed = track->PropagateInField(safety, kFALSE, trackin[0]);
      } 
      gNsafeSteps++;
      if (crossed) return crossed;  
      // Recompute snext and safety for todo tracks
      ComputeTransportLengthSingle(trackin);
   }
   return crossed;
}

//______________________________________________________________________________
void GeantVolumeBasket::TransportSingle()
{
// Transport all particles in this basket one by one (empty basket).
   if (!fNtracks) return;
   // Main loop
   Int_t itrack, nout;   
   Int_t *particles = gPartInd[0]->GetArray();
   Int_t *partnext  = gPartNext[0]->GetArray();
   memcpy(particles, fIndex, fNtracks*sizeof(Int_t));
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
//   Int_t tid = nav->GetThreadId();
//   gPushed[tid]->ResetAllBits();
   Bool_t transported = kFALSE;
   Bool_t crossed = kFALSE;
   Int_t generation = 0;
   TGeoBranchArray a;
   TGeoBranchArray b;
   Int_t n10 = fNtracks/10;
   for (itrack=0; itrack<fNtracks; itrack++) {
      if (n10) {
         if ((itrack%n10) == 0) Printf("%i percent", Int_t(100*itrack/fNtracks));
      }
      // Skip neutral particles for the moment (!!!)
      if (!gTracks[particles[itrack]]->charge) continue;
      transported = kFALSE;
      // Init navigator
      GetBranchArray(gTracks[particles[itrack]])->UpdateNavigator(nav);
      gVolume = gGeoManager->GetCurrentVolume();
      // Physics step      
      if (gUsePhysics) PhysicsSelect(1,&particles[itrack],0);
      while (!transported) {
         crossed = kFALSE;
         a.InitFromNavigator(nav);
         // Loop all tracks to generate physics/geometry steps
         generation++;
         // Geometry snext and safety
         ComputeTransportLengthSingle(&particles[itrack]);
         crossed = PropagateTrack(&particles[itrack]);
         if (!crossed) {
            // Do post-step actions
            if (gUsePhysics) {
               // Apply continuous processes first
               for (Int_t iproc=0; iproc<gNprocesses; iproc++) {
                  if (gProcessType[iproc] == kDiscrete) continue;
                  nout = 0;
                  PostStepAction[iproc](1, &particles[itrack], nout, partnext,0);
                  if (!nout) {
                     transported = kTRUE;
                     break;
                  }   
               }
               // Apply discrete process if selected
               if (!transported && gProcessType[gTracks[particles[itrack]]->process] == kDiscrete) {
                  nout = 0;
                  PostStepAction[gTracks[particles[itrack]]->process](1, &particles[itrack], nout, partnext,0);
                  if (!nout) {
                     transported = kTRUE;
                     break;
                  } 
               }     
               if (!transported) PhysicsSelect(1,&particles[itrack],0);
            }   
            continue;
         } else {
            // Energy deposition up to boundary
            if (gUsePhysics) {
               for (Int_t iproc=0; iproc<gNprocesses; iproc++) {
                  if (gProcessType[iproc] == kDiscrete) continue;
                  // In future take into account that particles may be generated 
                  // here (like dray)
                  Int_t nafter = 0;
                  PostStepAction[iproc](1, &particles[itrack], nafter, NULL,0);
                  ResetStep(1, &particles[itrack]);
               }   
            }
         }   
         // Exit setup ?
         if (nav->IsOutside()) break;
         // Particle crossed boundary
         b.InitFromNavigator(nav);
         if (b==a) continue; // Particle entered the same volume, just continue the step
         // Physics step      
         gVolume = gGeoManager->GetCurrentVolume();
         if (gUsePhysics) PhysicsSelect(1,&particles[itrack],0);
      }
   }
}   

//______________________________________________________________________________
void *TransportTracks(void *)
{
// Thread propagating all tracks from a basket.
   Int_t tid = TGeoManager::ThreadId();
//   Printf("(%d) WORKER started", tid);
   // Create navigator if none serving this thread.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   if (!nav) nav = gGeoManager->AddNavigator();

   Int_t indmin, indmax;
   Int_t ntotnext, ntmp, ntodo, ncross, cputime, ntotransport;
   GeantTrack *track = 0;
   Int_t itrack;
   Int_t *particles = 0;
   Int_t *partnext  = 0;
   Int_t *parttodo  = 0;
   Int_t *partcross = 0;
   Int_t generation = 0;
   GeantVolumeBasket *basket = 0;
   while (gCurrentBasket) {
      feeder_queue->wait_and_pop(ntmp);
//      Printf("Popped %d\n", ntmp);
//      bufferStart->Receive();
      basket = gCurrentBasket;
      ntotransport = basket->GetNtracks();  // all tracks to be transported      
      if (!ntotransport) goto finish;
      // Work splitting per thread
      basket->GetWorkload(indmin, indmax);
      ntotransport = indmax-indmin;
      if (!ntotransport) goto finish;      
//      Printf("(%d) ================= BASKET %s: %d tracks (%d-%d)", tid, basket->GetName(), ntotransport, indmin,indmax);
      particles = gPartInd[tid]->GetArray();
      partnext  = gPartNext[tid]->GetArray();
      parttodo  = gPartTodo[tid]->GetArray();
      partcross = gPartCross[tid]->GetArray();
      memcpy(particles, &basket->GetIndArray()[indmin], ntotransport*sizeof(Int_t));
//      PrintParticles(particles, ntotransport, tid);
      ntotnext = 0;
      ntmp = 0;
      ntodo = 0;
      ncross = 0;
      cputime = 0.;   
      generation = 0;
      track = 0;
      while (ntotransport) {
         generation++;
         // Loop all tracks to generate physics/geometry steps
         // Physics step
         if (gUsePhysics) PhysicsSelect(ntotransport, particles, tid);
         // Geometry snext and safety
         basket->ComputeTransportLength(ntotransport, particles);
         // Propagate tracks with physics step and check if they survive boundaries 
         // or physics
         ntmp = ntotransport;
         ntodo = ntotransport;
         ntotnext = 0;
         Int_t *ptrParticles = particles;
         // Propagate all tracks alive with physics step.
         // If a boundary is encountered before the physics step, stop the track
         if (gUseDebug) Printf("(%d) --- propagating %d tracks for volume basket %s", tid, ntodo, basket->GetName());
         while (ntodo) {
            ntodo = 0;
            ncross = 0;
            // Propagate ALL ntmp tracks
            basket->PropagateTracks(ntmp, ptrParticles, ntotnext, partnext, ntodo, parttodo, ncross, partcross);
//            printf("(%d) %s   ->crossing particles (%d): ",tid, basket->GetName(), ncross); 
//            for (Int_t ii=0; ii<ncross; ii++) printf("%d ", partcross[ii]);
//            printf("\n(%d) %s   ->remaining particles (%d): ", tid, basket->GetName(), ntotnext);
//            for (Int_t ii=0; ii<ntotnext; ii++) printf("%d ", partnext[ii]);
//            printf("\n(%d) %s   ->todo particles (%d): ", tid, basket->GetName(),ntodo);
//            for (Int_t ii=0; ii<ntodo; ii++) printf("%d ", parttodo[ii]);
//            printf("\n");
            // Post-step actions by continuous processes for particles reaching boundaries
            if (gUsePhysics && ncross) {
               for (Int_t iproc=0; iproc<gNprocesses; iproc++) {
                  if (gProcessType[iproc] == kDiscrete) continue;
                  Int_t nafter = 0;
                  PostStepAction[iproc](ncross, partcross, nafter, NULL,tid);
                  ResetStep(ncross, partcross);
               }   
            }      
            ntmp = ntodo;
            ptrParticles = parttodo;
         }
        
         // Copy only tracks that survived boundaries (well we will have to think of
         // those too, like passing them to the next volume...)
         memcpy(particles, partnext, ntotnext*sizeof(Int_t));
         ntotransport = ntotnext;
            
         // Do post-step actions on remaining particles
         ntotnext = 0;
         // Loop all processes to group particles per process
         if (gUsePhysics && ntotransport) {
            // Apply continuous processes to all particles
            for (Int_t iproc=0; iproc<gNprocesses; iproc++) {
               if (gProcessType[iproc] == kDiscrete) continue;
               ntodo = 0;
               PostStepAction[iproc](ntotransport, particles, ntodo, parttodo, tid);
               // Do we have stopped particles ?
               if (ntodo<ntotransport) {
                  memcpy(particles, parttodo, ntodo*sizeof(Int_t));
                  ntotransport = ntodo;
               }
            } 
            // Copy al tracks for which step was limited by a continuous process
            // to the next array
            for (Int_t itr=0; itr<ntotransport; itr++) {
               if (gProcessType[gTracks[particles[itr]]->process] == kContinuous)
                  partnext[ntotnext++] = particles[itr];
            }      
            // Discrete processes only
            for (Int_t iproc=0; iproc<gNprocesses; iproc++) {
               // Make arrays of particles per process -> ntodo, parttodo
               if (gProcessType[iproc] == kContinuous) continue;
               ntodo = 0;
               SelectTracksForProcess(iproc, ntotransport, particles, ntodo, parttodo);
               if (!ntodo) continue;
               if (gPartTodo[tid]->GetSize()-ntodo<500) {
                  gPartTodo[tid]->Set(2*gPartTodo[tid]->GetSize());
                  parttodo  = gPartTodo[tid]->GetArray(); 
               }   
               // Do post step actions for particles suffering a given process.
               // Surviving particles are added to the next array
      //         Printf("PostStep for proc %d: %d particles:\n", iproc, ntodo);
               PostStepAction[iproc](ntodo, parttodo, ntotnext, partnext,tid);
               if (gPartNext[tid]->GetSize()-ntotnext<500) {
                  gPartNext[tid]->Set(2*gPartNext[tid]->GetSize());
                  partnext  = gPartNext[tid]->GetArray();
                  gPartInd[tid]->Set(2*gPartInd[tid]->GetSize());
                  particles = gPartInd[tid]->GetArray();
               }   
            }
            memcpy(particles, partnext, ntotnext*sizeof(Int_t));
            ntotransport = ntotnext;
         }
         // I/O: Dump current generation
//         Printf("   ### Generation %d:  %d tracks  cputime=%f", generation, ntotransport,cputime);
         if (gFillTree) {
            cputime = gTimer.CpuTime();
            gOutdata->SetStamp(basket->GetVolume()->GetNumber(), gBasketGeneration, generation, ntotransport, cputime);
            for (itrack=0; itrack<ntotransport;itrack++) {
               track = gTracks[particles[itrack]];
               gOutdata->SetTrack(itrack, track);
            }   
            gOutTree->Fill();
         }
      }   

finish:
      // the last thread to finish wakes up the main thread
      // ... then go to sleep
      // Checkpoint. 
//      Printf("Thread %d finished", tid);
      answer_queue->push(tid);
//      bufferStop->StartN();
   }
   return 0;
}
   
//______________________________________________________________________________
void GeantVolumeBasket::Print(Option_t *) const
{
// Print info about the basket content.
   Printf("Volume basket %s: ntracks=%d npaths=%d", GetName(), GetNtracks(), fNpaths);
//   if (gUseDebug) {
      for (Int_t i=0; i<fNtracks; i++) {
//         if (gDebug && (gDebugTrk==fIndex[i] || gDebugTrk<0)) {
//            Printf("   %d - track %d: ", i, fIndex[i]);
            gTracks[fIndex[i]]->Print();
            if (gTracks[fIndex[i]]->path->GetNode(gTracks[fIndex[i]]->path->GetLevel())->GetVolume() != fVolume) {
               Printf("ERROR: Wrong path for basket %s", GetName());
               *((int*)0)=0;
            }  
//            GetBranchArray(i)->Print();
//         }   
      }
//   }      
}

//______________________________________________________________________________
void SortBaskets(Int_t *index)
{
// Sort baskets in decreasing number of particles. The order is set in the provided index array of size gNbaskets minimum.
   if (gNbaskets<1) return;
   Int_t *ipart = new Int_t[gNbaskets];
   ipart[0] = 0;
   for (Int_t ibasket=0; ibasket<gNbaskets; ibasket++) ipart[ibasket] = gBasketArray[ibasket]->GetNtracks();
   TMath::Sort(gNbaskets, ipart, index, kTRUE);
   delete [] ipart;
}   

// ######## Main loop

//______________________________________________________________________________
void PropagatorGeom(const char *geomfile="geometry.root", Int_t nthreads=4, Bool_t graphics=kFALSE, Bool_t single=kFALSE, Double_t vertx=0., Double_t verty=0., Double_t vertz=0.)
{
// Propagate gNevents in the volume containing the vertex. 
// Simulate 2 physics processes up to exiting the current volume.
   static Bool_t called=kFALSE;
   Int_t ipop;
   gNthreads = nthreads;
   gSingleTrack = single;
   if (called) {
      Printf("Sorry, you can call this only once per session.");
      return;
   }
   called = kTRUE;   
   // Initialize pointers to functions
   InitFunctions();
   gKineTF1 = new TF1("gKineTF1","gaus",emin,emax);
   gKineTF1->SetParameters(1,3*emin,5);
   // Initialize vertex
   gVertex[0] = vertx;
   gVertex[1] = verty;
   gVertex[2] = vertz;
//   Int_t itrack;

   // Initialize geometry and current volume
   if (!LoadGeometry(geomfile)) return;
   if (gSingleTrack) Printf("==== Executing in single track loop mode using %d threads ====", gNthreads);
   else              Printf("==== Executing in vectorized mode using %d threads ====",gNthreads);
   if (gFillTree)    Printf("  I/O enabled - disable if comparing single track loop with vectorized modes");
   else              Printf("  I/O disabled");
   if (gUsePhysics)  Printf("  Physics ON with %d processes", gNprocesses);
   else              Printf("  Physics OFF");
   // Create the basket array
   Int_t nvols = gGeoManager->GetListOfVolumes()->GetEntries();
   gBasketArray = new GeantVolumeBasket*[nvols];

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
   GeantVolumeBasket *basket = ImportTracks(gNevents, gNaverage);
   gBasketArray[gNbaskets++] = basket;

   // Initialize tree
   gOutdata = new GeantOutput();
   gOutdata->Init(gMaxTracks);
   gOutFile = 0;
   if (gFillTree) {
      gOutFile = new TFile("output.root", "RECREATE");
      gOutTree = new TTree("TK","Transport track data");
      gOutTree->Branch("gen", &gOutdata);
   }
   
   // Initialize threads
   TThread *t;
   TList *listThreads = new TList();
   listThreads->SetOwner();
   for (Int_t ith=0; ith<gNthreads; ith++) {
      t = new TThread(TransportTracks);
      listThreads->Add(t);
   } 
   Bool_t threadsStarted = kFALSE;
   
   
   // Loop baskets and transport particles until there is nothing to transport anymore
   gTransportOngoing = kTRUE;
   gGeoManager->SetMultiThread(kTRUE);
   Int_t nbaskets, nb0, ntrackgen;
   gBasketGeneration = 0;
   gTimer.Start();
   while (gTransportOngoing) {
      if (gSingleTrack) {
         basket->TransportSingle();
         break;
      }    
      gTransportOngoing = kFALSE;
      nbaskets = gNbaskets;
      // Loop current generation of baskets
      Int_t *index = new Int_t[nbaskets];
      SortBaskets(index);
      Bool_t useThreshold = kFALSE;
      Int_t nbtrue = 0;
      if (gBasketArray[index[0]]->GetNtracks()>gNminThreshold) useThreshold = kTRUE;
      nb0 = 0;
      ntrackgen = 0;
      for (Int_t ibasket=0; ibasket<nbaskets; ibasket++) {
         Int_t ntracks = gBasketArray[index[ibasket]]->GetNtracks();
         if (!ntracks) continue;
         if (useThreshold && ntracks<gNminThreshold) continue;
         ntrackgen += ntracks;
         nb0++;
      }   
      if (nbaskets) Printf("#### GENERATION #%04d (%05d part) OF %04d/%04d VOLUME BASKETS, (TOP= %s) ####", gBasketGeneration, ntrackgen, nb0, nbaskets, gBasketArray[index[0]]->GetName());
      for (Int_t ibasket=0; ibasket<nbaskets; ibasket++) {
         gCurrentBasket = gBasketArray[index[ibasket]];
         Int_t ntracks = gCurrentBasket->GetNtracks();
         if (!ntracks) continue;
         if (useThreshold && ntracks<gNminThreshold) {
            if (ntracks) gTransportOngoing=kTRUE;
            continue;
         }   
         if (graphics) {
            if (ibasket < 102) {
               hbaskets->Fill(ibasket,ntracks);
               hbaskets->SetTitle(Form("baskets population for generation %d, volume = %s",gBasketGeneration,gCurrentBasket->GetName()));
               pad2->Modified();
               c1->Update();
            }
         } else {
//            gCurrentBasket->Print();
         }
         nbtrue++;
         // Start threaded transport
//         gCurrentBasket->TransportTracks();
//         Printf("CURRENT BASKET: %s", gCurrentBasket->GetName());
         gVolume = gCurrentBasket->GetVolume();
//         gCurrentBasket->Print();
         if (!threadsStarted) {
            for (Int_t ith=0; ith<gNthreads; ith++) {
               t = (TThread*)listThreads->At(ith);
               t->Run();
            }
//         TThread::Sleep(0,2000000);
            threadsStarted = kTRUE;
         }
         // Put the object in buffer for N threads
//         bufferStart->Start();
         Int_t nchunk = ntracks/gNthreads;
         Int_t nworkers = gNthreads;
         if (!nchunk) nworkers = ntracks;
         for (Int_t iwork=0; iwork<nworkers; iwork++) feeder_queue->push(iwork);
//         Printf("== %d objects put by main thread", nworkers);
         // Retreive the result
         while(nworkers) {
            answer_queue->wait_and_pop(ipop);
//            Printf("Worker %d finished", ipop);
            nworkers--;
         }   
//         bufferStop->ReceiveN();
         gCurrentBasket->Clear();
//         Printf("== basket cleared");
//         TThread::Sleep(0,2000000);
      }
      delete [] index;
      if (graphics) {
         hnb->Fill(gBasketGeneration,nbtrue);
         pad1->Modified();
         hbaskets->Reset();
      }
      gBasketGeneration++;
   }      
   gTimer.Stop();
   gTimeCounter.Print();
//   for (Int_t itr=0; itr<gNtracks; itr++) gTracks[itr]->Print();
   if (gFillTree) gOutTree->AutoSave();
   delete gOutFile;
   Double_t rtime = gTimer.RealTime();
   Double_t ctime = gTimer.CpuTime();
   gTimer.Print();
   const char *geomname=geomfile;
   if(strstr(geomfile,"http://root.cern.ch/files/")) geomname=geomfile+strlen("http://root.cern.ch/files/");
   Printf("=== Transported: %lld,  safety steps: %lld,  snext steps: %lld, RT=%gs, CP=%gs", gNtransported, gNsafeSteps, gNsnextSteps,rtime,ctime);
   gSystem->mkdir("results");
   FILE *fp = fopen(Form("results/%s_%d.dat",geomname,single),"w");
   fprintf(fp,"%d %lld %lld %g %g",single, gNsafeSteps, gNsnextSteps,rtime,ctime);
   fclose(fp);
   delete listThreads;
   gOutFile = 0;
   gOutTree = 0;
}

//______________________________________________________________________________
Bool_t DrawData(Int_t color = kRed)
{
// Draw a given generation of track points.
   static Long64_t ientry = 0;
   Long64_t nentries = gOutTree->GetEntries();
   if (ientry==nentries) return kFALSE;   
   TPolyMarker3D *pmgen = new TPolyMarker3D();
   pmgen->SetMarkerColor(color);
   Int_t nread = 0;
   Int_t ibgen = 0;
   while (ientry<nentries) {
      gOutTree->GetEntry(ientry++);
      if (!nread++) {
         ibgen = gOutdata->fBasketGeneration;
      }   
      if (gOutdata->fBasketGeneration > ibgen) {
         ientry--;
         break;
      }      
      for (Int_t itrack=0; itrack<gOutdata->fNtracks; itrack++) 
         pmgen->SetNextPoint(gOutdata->fX[itrack], gOutdata->fY[itrack],gOutdata->fZ[itrack]);
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
   if (!gOutFile) gOutFile = new TFile("output.root");
   if (!gOutTree) {
      gOutTree = (TTree*)gOutFile->Get("TK");
      gOutdata = new GeantOutput();
      gOutTree->SetBranchAddress("gen", &gOutdata);
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
   if (!gGeoManager) LoadGeometry(file);
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
