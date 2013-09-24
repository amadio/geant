#ifndef GEANT_TRACK
#define GEANT_TRACK

#include "globals.h"
#include "TMath.h"
#include "TBits.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h> 
#else
#include "mm_malloc.h"
#endif

#ifndef ALIGN_PADDING
#define ALIGN_PADDING 32 
#endif

class TGeoBranchArray;
class GeantMainScheduler;

const Double_t kB2C = -0.299792458e-3;
enum TrackStatus_t {kAlive, kKilled, kBoundary, kExitingSetup, kPhysics, kPostponed};
enum TransportAction_t {
   kPostpone = 0,   // return imediately and postpone whatever tracks left
   kSingle   = 1,   // perform remaining loop in single track mode
   kVector   = 2    // perform remaining loop in vectorized mode
};   

// Rounding up the desired allocation value to the closest padding multiple
int round_up_align(int num)
{
   int remainder = num % ALIGN_PADDING;
   if (remainder == 0) return num;
   return (num+ALIGN_PADDING-remainder);
}  

class GeantTrack_v;

//______________________________________________________________________________
class GeantTrack : public TObject {
private:
   Int_t    fEvent;     // event number
   Int_t    fEvslot;    // event slot
   Int_t    fParticle;  // index of corresponding particle
   Int_t    fPDG;       // particle pdg code
   Int_t    fG5code;    // G5 particle code
   Int_t    fCharge;    // particle charge
   Int_t    fProcess;   // current process
   Int_t    fIzero;     // number of small steps used to catch errors
   Int_t    fNsteps;    // number of steps made
   Species_t fSpecies;  // particle species
   TrackStatus_t fStatus; // track status
   Double_t fMass;      // particle mass
   Double_t fXpos;      // position
   Double_t fYpos;
   Double_t fZpos;
   Double_t fXdir;      // direction
   Double_t fYdir;
   Double_t fZdir;
   Double_t fP;         // momentum
   Double_t fE;         // energy
   Double_t fPstep;     // selected physical step
   Double_t fStep;      // current step
   Double_t fSnext;     // straight distance to next boundary
   Double_t fSafety;    // safe distance to any boundary
   Bool_t   fFrombdr;   // true if starting from boundary
   Bool_t   fPending;
   TGeoBranchArray *fPath; // path for this particle in the geometry
   TGeoBranchArray *fNextpath; // path for next volume
   
public:
   GeantTrack();
   GeantTrack(const GeantTrack &other);
   GeantTrack &operator=(const GeantTrack &other);
   GeantTrack(Int_t ipdg);
   ~GeantTrack();

   Double_t           Beta()  const {return fP/fE;}
   Int_t              Charge() const {return fCharge;}
   Double_t           Curvature() const {return (fCharge)?TMath::Abs(kB2C*gPropagator->fBmag/Pt()):0.;}
   const Double_t    *Direction() const {return &fXdir;}
   Double_t           DirX() const {return fXdir;}
   Double_t           DirY() const {return fYdir;}
   Double_t           DirZ() const {return fZdir;}
   Double_t           E() const {return fE;}
   Int_t              Event() const {return fEvent;}
   Int_t              EventSlot() const  {return fEvslot;}
   Bool_t             FromBoundary() const {return fFrombdr;}
   Int_t              G5code() const {return fG5code;}
   Double_t           Gamma() const {return fMass?fE/fMass:TMath::Limits<double>::Max();}
   Double_t           GetPstep() const {return fPstep;}
   TGeoBranchArray   *GetPath() const {return fPath;}
   TGeoBranchArray   *GetNextPath() const {return fNextpath;}
   Int_t              GetNsteps() const {return fNsteps;}
   Double_t           GetStep() const {return fStep;}
   Double_t           GetSnext() const {return fSnext;}
   Double_t           GetSafety() const {return fSafety;}
   Bool_t             IsAlive() const {return (fStatus != kKilled);}
   Bool_t             IsOnBoundary() const {return (fStatus == kBoundary);}
   Int_t              Izero() const {return fIzero;}
   void               Kill()        {fStatus = kKilled;}
   Double_t           Mass() const {return fMass;}
   Double_t           P() const {return fP;}
   Double_t           Px() const {return fP*fXdir;}
   Double_t           Py() const {return fP*fYdir;}
   Double_t           Pz() const {return fP*fZdir;}
   Double_t           Pt()    const {return fP*TMath::Sqrt(fXdir*fXdir+fYdir*fYdir);}
   Int_t              Particle() const {return fParticle;}
   Bool_t             Pending() const {return fPending;}
   Int_t              PDG() const   {return fPDG;}
   Int_t              Process() const {return fProcess;}
   const Double_t    *Position() const {return &fXpos;}
   Double_t           PosX() const {return fXpos;}
   Double_t           PosY() const {return fYpos;}
   Double_t           PosZ() const {return fZpos;}
   void               Print(Int_t trackindex=0) const;
   Bool_t             PropagateInFieldSingle(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateInField(Double_t step, Bool_t checkcross, Int_t itr);
   GeantVolumeBasket *PropagateStraight(Double_t step, Int_t itrack);
   Species_t          Species() const {return fSpecies;}
   TrackStatus_t      Status() const  {return fStatus;}
   
   void               Reset();
   Double_t           X() const {return fXpos;}
   Double_t           Y() const {return fYpos;}
   Double_t           Z() const {return fZpos;}
   Bool_t             ZeroCounter() const {return fIzero;}
   
   void               ReadFromVector(const GeantTrack_v &arr, Int_t i);
   void               SetEvent(Int_t event) {fEvent = event;}
   void               SetEvslot(Int_t slot) {fEvslot = slot;}
   void               SetParticle(Int_t particle) {fParticle = particle;}
   void               SetPDG(Int_t pdg) {fPDG = pdg;}
   void               SetG5code(Int_t g5code) {fG5code = g5code;}
   void               SetCharge(Int_t charge) {fCharge = charge;} 
   void               SetProcess(Int_t process) {fProcess = process;}
   void               SetIzero(Int_t izero) {fIzero = izero;}
   void               SetNsteps(Int_t nsteps) {fNsteps = nsteps;}
        
   
   ClassDef(GeantTrack, 1)      // The track
};

//______________________________________________________________________________
// SOA for GeantTrack used at processing time

struct GeantTrack_v {
public:
   Int_t     fNtracks;    // number of tracks contained
   Int_t     fMaxtracks;  // max size for tracks
   Int_t     fNselected;  // Number of selected tracks
   TBits     fHoles;      // Bits of holes
   TBits     fSelected;   // Mask of selected tracks for the current operation
   Bool_t    fCompact;    // Flag marking the compactness
   char     *fBuf;        // buffer holding tracks data

   Int_t    *fEventV;     // event numbers
   Int_t    *fEvslotV;    // event slots
   Int_t    *fParticleV;  // indices of corresponding particles
   Int_t    *fPDGV;       // particle pdg codes
   Int_t    *fG5codeV;    // G5 internal codes
   Int_t    *fChargeV;    // particle charges
   Int_t    *fProcessV;   // current process
   Int_t    *fIzeroV;     // number of small steps used to catch errors
   Int_t    *fNstepsV;    // number of steps made
   Species_t *fSpeciesV;  // particle species
   TrackStatus_t *fStatusV; // track statuses
   Double_t *fMassV;      // particle masses
   Double_t *fXposV;    // Array of track positions
   Double_t *fYposV;
   Double_t *fZposV;
   Double_t *fXdirV;    // Array of track directions
   Double_t *fYdirV;
   Double_t *fZdirV;
   Double_t *fPV;         // momenta
   Double_t *fEV;         // energies
   Double_t *fPstepV;     // selected physical steps
   Double_t *fStepV;      // current steps
   Double_t *fSnextV;     // straight distances to next boundary
   Double_t *fSafetyV;     // safe distances to any boundary
   Bool_t   *fFrombdrV;     // true if starting from boundary
   Bool_t   *fPendingV;
   TGeoBranchArray **fPathV; // paths for the particles in the geometry
   TGeoBranchArray **fNextpathV; // paths for next volumes

   void AssignInBuffer(const char *buff, Int_t size);
   void CopyToBuffer(const char *buff, Int_t size);

public:   
   GeantTrack_v();
   GeantTrack_v(Int_t size);
   GeantTrack_v(const GeantTrack_v &track_v);
   GeantTrack_v &operator=(const GeantTrack_v &track_v);
   virtual ~GeantTrack_v();

   Int_t     Capacity() const     {return fMaxtracks;}
   Int_t     GetNtracks() const   {return fNtracks;}
   Int_t     GetNselected() const {return fNselected;}
   void      AddTrack(const GeantTrack &track);
   void      AddTrack(const GeantTrack_v &arr, Int_t i);
   void      AddTracks(const GeantTrack_v &arr, Int_t istart, Int_t iend);
   Int_t     FlushTracks(GeantMainScheduler *main);
   void      MarkRemoved(Int_t i) {fHoles.SetBitNumber(i); fCompact=kFALSE;}
   void      RemoveTracks(Int_t from, Int_t to);
   void      Deselect(Int_t i)    {fSelected.SetBitNumber(i, kFALSE);}
   void      Select(Int_t i)      {fSelected.SetBitNumber(i);}
   void      SelectTracks(Int_t n) {fNselected = n;}
   Bool_t    IsSelected(Int_t i)  {return fSelected.TestBitNumber(i);}
   virtual void      Clear(Option_t *option="");
   Int_t     Compact(GeantTrack_v *moveto=0);
   Bool_t    Contains(Int_t evstart, Int_t nevents=1) const;
   void      ClearSelection()     {fSelected.ResetAllBits();}
   void      GetTrack(Int_t i, GeantTrack &track) const;
   Bool_t    IsCompact() const {return fCompact;}
      
   void Print() {
      printf("fXposV=%p fYposV=%p fZposV=%p fXdirV=%p fYdirV=%p fZdirV=%p ptot=%p fG5codeV=%p\n",
              fXposV,fYposV,fZposV,fXdirV,fYdirV,fZdirV,fPV,fG5codeV);
   }
   void      NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep, 
                       const Double_t *x, const Double_t *y, const Double_t *z,
                       const Double_t *dirx, const Double_t *diry, const Double_t *dirz,
                       TGeoBranchArray **pathin, TGeoBranchArray **pathout, 
                       Double_t *step, Double_t *safe, Bool_t *isonbdr);
   void      NavIsSameLocation(Int_t ntracks, TGeoBranchArray **start, TGeoBranchArray **end, Bool_t *same);
   Bool_t    NavIsSameLocationSingle(Int_t itr, TGeoBranchArray **start, TGeoBranchArray **end);
   TransportAction_t PostponedAction() const;
   Int_t     PostponeTracks(GeantTrack_v &output);
   void      PropagateBack(Int_t itr, Double_t crtstep);
   Int_t     PropagateInField(Int_t ntracks, const Double_t *crtstep);
   Int_t     PropagateInFieldSingle(Int_t itr, Double_t crtstep, Bool_t checkcross);
   void      ComputeTransportLength(Int_t ntracks);
   void      ComputeTransportLengthSingle(Int_t itr);
   void      PropagateInVolume(const Double_t *crtstep);
   void      PropagateInVolumeSingle(Int_t i, Double_t crtstep);
   Int_t     PropagateStraight(Int_t ntracks, Double_t *crtstep);
   Int_t     PropagateTracks(GeantTrack_v &output);
   Int_t     PropagateTracksSingle(GeantTrack_v &output, Int_t stage);
   
   void      Resize(Int_t newsize);
   void      ReplaceTrack(Int_t i, Int_t withj);
   Int_t     Reshuffle();
   void      SwapTracks(Int_t i, Int_t j);
// Track methods
   Double_t           Beta(Int_t i)  const {return fPV[i]/fEV[i];}
   Double_t           Curvature(Int_t i) const {return (fChargeV[i])?TMath::Abs(kB2C*gPropagator->fBmag/Pt(i)):0.;}
   Double_t           Gamma(Int_t i) const {return fMassV[i]?fEV[i]/fMassV[i]:TMath::Limits<double>::Max();}
   Double_t           Px(Int_t i) const {return fPV[i]*fXdirV[i];}
   Double_t           Py(Int_t i) const {return fPV[i]*fYdirV[i];}
   Double_t           Pz(Int_t i) const {return fPV[i]*fZdirV[i];}
   Double_t           Pt(Int_t i) const {return fPV[i]*TMath::Sqrt(fXdirV[i]*fXdirV[i]+fYdirV[i]*fYdirV[i]);}

   ClassDef(GeantTrack_v, 1)      // SOA for GeantTrack class           
};

#endif
