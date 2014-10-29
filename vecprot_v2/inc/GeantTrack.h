#ifndef GEANT_TRACK
#define GEANT_TRACK

//#include "globals.h"
#include "TMath.h"
#include "TBits.h"

#ifdef __STAT_DEBUG
#include "GeantTrackStat.h"
#endif   

#if __cplusplus >= 201103L
#include <atomic>
#endif

#ifndef ALIGN_PADDING
#define ALIGN_PADDING 32 
#endif

#ifdef USE_VECGEOM_NAVIGATOR
 #include "navigation/NavigationState.h"
 typedef VECGEOM_NAMESPACE::NavigationState VolumePath_t;
#else
 #include "TGeoBranchArray.h"       // needed due to templated pools
 typedef TGeoBranchArray VolumePath_t;
#endif


const Double_t kB2C = -0.299792458e-3;
enum TrackStatus_t {kAlive, kKilled, kInFlight, kBoundary, kExitingSetup, kPhysics, kPostponed, kNew};
enum TransportAction_t {
   kDone     = 0,   // Return immediately - no tracks left
   kPostpone = 1,   // return imediately and postpone whatever tracks left
   kSingle   = 2,   // perform remaining loop in single track mode
   kVector   = 3    // perform remaining loop in vectorized mode
};   
// types
enum Species_t {kHadron, kLepton};

class GeantTrack_v;

//______________________________________________________________________________
class GeantTrack : public TObject {
public:
   Int_t    fEvent;     // event number
   Int_t    fEvslot;    // event slot
   Int_t    fParticle;  // index of corresponding particle
   Int_t    fPDG;       // particle pdg code
   Int_t    fG5code;    // G5 particle code
   Int_t    fEindex;    // element index
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
   Double_t fEdep;      // energy deposition in the step
   Double_t fPstep;     // selected physical step
   Double_t fStep;      // current step
   Double_t fSnext;     // straight distance to next boundary
   Double_t fSafety;    // safe distance to any boundary
   Bool_t   fFrombdr;   // true if starting from boundary
   Bool_t   fPending;
   VolumePath_t *fPath;
   VolumePath_t *fNextpath;
   
public:
   GeantTrack();
   GeantTrack(const GeantTrack &other);
   GeantTrack &operator=(const GeantTrack &other);
   GeantTrack(Int_t ipdg);
   ~GeantTrack();

   Double_t           Beta()  const {return fP/fE;}
   Int_t              Charge() const {return fCharge;}
   Double_t           Curvature() const;
   const Double_t    *Direction() const {return &fXdir;}
   Double_t           DirX() const {return fXdir;}
   Double_t           DirY() const {return fYdir;}
   Double_t           DirZ() const {return fZdir;}
   Double_t           E() const {return fE;}
   Double_t           Edep() const {return fEdep;}
   Int_t              Event() const {return fEvent;}
   Int_t              EventSlot() const  {return fEvslot;}
   Bool_t             FromBoundary() const {return fFrombdr;}
   Int_t              G5code() const {return fG5code;}
   Int_t              EIndex() const {return fEindex;}
   Double_t           Gamma() const {return fMass?fE/fMass:TMath::Limits<double>::Max();}
   Double_t           GetPstep() const {return fPstep;}
   VolumePath_t*      GetPath() const {return fPath;}
   VolumePath_t*      GetNextPath() const {return fNextpath;}
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
   Species_t          Species() const {return fSpecies;}
   TrackStatus_t      Status() const  {return fStatus;}
   
   virtual void       Clear(Option_t *option="");
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
   void               SetEindex(Int_t ind) {fEindex = ind;}
   void               SetCharge(Int_t charge) {fCharge = charge;} 
   void               SetProcess(Int_t process) {fProcess = process;}
   void               SetIzero(Int_t izero) {fIzero = izero;}
   void               SetNsteps(Int_t nsteps) {fNsteps = nsteps;}
   void               SetSpecies(Species_t   species) {fSpecies = species;}
   void               SetStatus(TrackStatus_t &status) {fStatus = status;}
   void               SetMass(Double_t mass) {fMass = mass;}
   void               SetPosition(Double_t x, Double_t y, Double_t z) {fXpos=x; fYpos=y; fZpos=z;}
   void               SetDirection(Double_t dx, Double_t dy, Double_t dz) {fXdir=dx; fYdir=dy; fZdir=dz;}
   void               SetP(Double_t p) {fP=p;}
   void               SetE(Double_t e) {fE = e;}
   void               SetEdep(Double_t edep) {fEdep = edep;}
   void               SetPstep(Double_t pstep) {fPstep = pstep;}
   void               SetSnext(Double_t snext) {fSnext = snext;}
   void               SetSafety(Double_t safety) {fSafety = safety;}
   void               SetFrombdr(Bool_t flag) {fFrombdr = flag;}
   void               SetPending(Bool_t flag) {fPending = flag;}
   void               SetPath(VolumePath_t const * const path);
   void               SetNextPath(VolumePath_t const * const path);
   
   ClassDef(GeantTrack, 1)      // The track
};

//______________________________________________________________________________
// SOA for GeantTrack used at processing time

class GeantTrack_v {
public:
#if __cplusplus >= 201103L
   std::atomic_int   fNtracks;  // number of tracks contained
#else
   Int_t     fNtracks;    // number of tracks contained  
#endif
   Int_t     fMaxtracks;  // max size for tracks
   Int_t     fNselected;  // Number of selected tracks
   TBits     fHoles;      // Bits of holes
   TBits     fSelected;   // Mask of selected tracks for the current operation
   Bool_t    fCompact;    // Flag marking the compactness
   
#ifdef __STAT_DEBUG_TRK
   GeantTrackStat fStat;  //! Statistics for the track container
#endif   
   Int_t     fMaxDepth;   // Maximum geometry depth allowed
   size_t    fBufSize;    // Size of the internal buffer
   char     *fVPstart;    // address of volume path buffer
   char     *fBuf;        // buffer holding tracks data

   Int_t    *fEventV;     // event numbers
   Int_t    *fEvslotV;    // event slots
   Int_t    *fParticleV;  // indices of corresponding particles
   Int_t    *fPDGV;       // particle pdg codes
   Int_t    *fG5codeV;    // G5 internal codes
   Int_t    *fEindexV;    // Element indices
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
   Double_t *fEdepV;      // Energy depositions
   Double_t *fPstepV;     // selected physical steps
   Double_t *fStepV;      // current steps
   Double_t *fSnextV;     // straight distances to next boundary
   Double_t *fSafetyV;     // safe distances to any boundary
   Bool_t   *fFrombdrV;     // true if starting from boundary
   Bool_t   *fPendingV;
   VolumePath_t **fPathV; // paths for the particles in the geometry
   VolumePath_t **fNextpathV; // paths for next volumes

   void AssignInBuffer(char *buff, Int_t size);
   void CopyToBuffer(char *buff, Int_t size);

private:
   GeantTrack_v(const GeantTrack_v &track_v); // not allowed
   GeantTrack_v &operator=(const GeantTrack_v &track_v); // not allowed

public:   
   GeantTrack_v();
   GeantTrack_v(Int_t size, Int_t maxdepth);
   virtual ~GeantTrack_v();

   size_t    BufferSize() const   {return fBufSize;}
   Int_t     Capacity() const     {return fMaxtracks;}
   static Bool_t IsSame(const GeantTrack_v &tr1, Int_t i1, const GeantTrack_v &tr2, Int_t i2);
#if __cplusplus >= 201103L
   Int_t     GetNtracks() const   {return fNtracks.load();}
   void      SetNtracks(Int_t ntracks) {fNtracks.store(ntracks);}
#else
   Int_t     GetNtracks() const   {return fNtracks;}
   void      SetNtracks(Int_t ntracks) {fNtracks = ntracks;}
#endif
   Int_t     GetNselected() const {return fNselected;}
#ifdef __STAT_DEBUG_TRK
   GeantTrackStat      &GetTrackStat() {return fStat;}
#endif   
   Int_t     AddTrack(GeantTrack &track, Bool_t import=kFALSE);
   Int_t     AddTrackSync(GeantTrack &track);
   Int_t     AddTrack(GeantTrack_v &arr, Int_t i, Bool_t import=kFALSE);
   Int_t     AddTrackSync(GeantTrack_v &arr, Int_t i);
   void      AddTracks(GeantTrack_v &arr, Int_t istart, Int_t iend, Bool_t import=kFALSE);
   void      CheckTracks();
   void      MarkRemoved(Int_t i) {fHoles.SetBitNumber(i); fCompact=kFALSE;}
   void      RemoveTracks(Int_t from, Int_t to);
   void      DeleteTrack(Int_t itr);
   void      Deselect(Int_t i)    {fSelected.SetBitNumber(i, kFALSE);}
   void      DeselectAll()        {fSelected.ResetAllBits(); fNselected = 0;}
   void      Select(Int_t i)      {fSelected.SetBitNumber(i);}
   void      SelectTracks(Int_t n) {fNselected = n;}
   Int_t     SortByStatus(TrackStatus_t status);
   Int_t     RemoveByStatus(TrackStatus_t status, GeantTrack_v &output);
   Bool_t    IsSelected(Int_t i)  {return fSelected.TestBitNumber(i);}
   virtual void      Clear(Option_t *option="");
   Int_t     Compact(GeantTrack_v *moveto=0);
   Bool_t    Contains(Int_t evstart, Int_t nevents=1) const;
   void      ClearSelection()     {fSelected.ResetAllBits();}
   void      GetTrack(Int_t i, GeantTrack &track) const;
   Bool_t    IsCompact() const {return fCompact;}
   void PrintPointers() {
      printf("fEventV=%p fFrombdrV=%p\n",  (void*)fEventV,(void*)fFrombdrV);
   }
   void PrintTrack(Int_t itr) const;
   void PrintTracks() const;

   void      NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep, 
                       const Double_t *x, const Double_t *y, const Double_t *z,
                       const Double_t *dirx, const Double_t *diry, const Double_t *dirz,
                       VolumePath_t **pathin, VolumePath_t **pathout,
                       Double_t *step, Double_t *safe, Bool_t *isonbdr, const GeantTrack_v *trk);
   void      NavIsSameLocation(Int_t ntracks, VolumePath_t **start, VolumePath_t **end, Bool_t *same);
   Bool_t    NavIsSameLocationSingle(Int_t itr, VolumePath_t **start, VolumePath_t **end);

   //void InspectGeometryState(Int_t itr) const;
   //void InspectIsSameLocation(Int_t itr) const;
#ifdef USE_VECGEOM_NAVIGATOR
    void CheckLocationPathConsistency(Int_t itr) const;
#endif

   TransportAction_t PostponedAction(Int_t ntracks) const;
   Int_t     PostponeTrack(Int_t itr, GeantTrack_v &output);
   Int_t     PostponeTracks(GeantTrack_v &output);
   //void      PropagateBack(Int_t itr, Double_t crtstep);
   void      ComputeTransportLength(Int_t ntracks);
   void      ComputeTransportLengthSingle(Int_t itr);
   void      PropagateInVolume(Int_t ntracks, const Double_t *crtstep, Int_t tid);
   void      PropagateInVolumeSingle(Int_t i, Double_t crtstep, Int_t tid);
   Int_t     PropagateStraight(Int_t ntracks, Double_t *crtstep);
   Int_t     PropagateTracks(GeantTrack_v &output, Int_t tid);
   Int_t     PropagateTracksSingle(GeantTrack_v &output, Int_t tid, Int_t stage=0);
   
   void      Resize(Int_t newsize);
   void      ReplaceTrack(Int_t i, Int_t withj);
   Int_t     Reshuffle();
   void      SwapTracks(Int_t i, Int_t j);
// Track methods
   Double_t           Beta(Int_t i)  const {return fPV[i]/fEV[i];}
   Double_t           Curvature(Int_t i) const;
   Double_t           SafeLength(Int_t i, Double_t eps=1.E-4);
   Double_t           Gamma(Int_t i) const {return fMassV[i]?fEV[i]/fMassV[i]:TMath::Limits<double>::Max();}
   Double_t           Px(Int_t i) const {return fPV[i]*fXdirV[i];}
   Double_t           Py(Int_t i) const {return fPV[i]*fYdirV[i];}
   Double_t           Pz(Int_t i) const {return fPV[i]*fZdirV[i];}
   Double_t           Pt(Int_t i) const {return fPV[i]*TMath::Sqrt(fXdirV[i]*fXdirV[i]+fYdirV[i]*fYdirV[i]);}
   static Int_t round_up_align(Int_t num) {
      int remainder = num % ALIGN_PADDING;
      if (remainder == 0) return num;
      return (num+ALIGN_PADDING-remainder);
   }  

   ClassDef(GeantTrack_v, 1)      // SOA for GeantTrack class           
};

#endif
