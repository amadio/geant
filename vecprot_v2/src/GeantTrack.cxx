#include "globals.h"
#include "TGeoBranchArray.h"
#include "TGeoNavigator.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoHelix.h"
#include "GeantTrack.h"
#include "GeantVolumeBasket.h"
#include "GeantThreadData.h"
#include "WorkloadManager.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h> 
#else
#include "mm_malloc.h"
#endif


const Double_t gTolerance = TGeoShape::Tolerance();

ClassImp(GeantTrack)

//______________________________________________________________________________
GeantTrack::GeantTrack() 
           :TObject(),
            fEvent(-1),
            fEvslot(-1),
            fParticle(-1),
            fPDG(0),
            fG5code(0),
            fCharge(0),
            fProcess(-1),
            fIzero(0), 
            fNsteps(0),
            fSpecies(kHadron),
            fStatus(kAlive),
            fMass(0),
            fXpos(0),
            fYpos(0),
            fZpos(0),
            fXdir(0),
            fYdir(0),
            fZdir(0),
            fP(0),
            fE(0),
            fPstep(1.E20), 
            fStep(0), 
            fSnext(0), 
            fSafety(0), 
            fFrombdr(false), 
            fPending(false)
            fPath(0),
            fNextpath(0),
{
// Dummy constructor
}

//______________________________________________________________________________
GeantTrack::GeantTrack(Int_t ipdg) 
           :TObject(),
            fEvent(-1),
            fEvslot(-1),
            fParticle(-1),
            fPDG(ipdg),
            fG5code(0),
            fCharge(0),
            fProcess(-1),
            fIzero(0), 
            fNsteps(0),
            fSpecies(kHadron),
            fStatus(kAlive),
            fMass(0),
            fXpos(0),
            fYpos(0),
            fZpos(0),
            fXdir(0),
            fYdir(0),
            fZdir(0),
            fP(0),
            fE(0),
            fPstep(1.E20), 
            fStep(0), 
            fSnext(0), 
            fSafety(0), 
            fFrombdr(false), 
            fPending(false),
            fPath(0),
            fNextpath(0)
{
// Constructor
   fPath = new TGeoBranchArray(30);
   fNextpath = new TGeoBranchArray(30);
}

//______________________________________________________________________________
GeantTrack::GeantTrack(const GeantTrack& other)
           :TObject(other),
            fEvent(other.fEvent),
            fEvslot(other.fEvslot),
            fParticle(other.fParticle),
            fPDG(other.fPDG),
            fG5code(other.fG5code),
            fCharge(other.fCharge),
            fProcess(other.fProcess),
            fIzero(other.fIzero), 
            fNsteps(other.fNsteps),
            fSpecies(other.fSpecies),
            fStatus(other.fStatus),
            fMass(other.fMass),
            fXpos(other.fXpos),
            fYpos(other.fYpos),
            fZpos(other.fZpos),
            fXdir(other.fXdir),
            fYdir(other.fYdir),
            fZdir(other.fZdir),
            fP(other.fP),
            fE(other.fE),
            fPstep(other.fPstep), 
            fStep(other.fStep), 
            fSnext(other.fSnext), 
            fSafety(other.fSafety), 
            fFrombdr(other.fFrombdr), 
            fPending(other.fPending),
            fPath(new TGeoBranchArray(*other.fPath)),
            fNextpath(new TGeoBranchArray(*other.fNextpath))
{
// Copy constructor
}

//______________________________________________________________________________
GeantTrack & GeantTrack::operator=(const GeantTrack &other)
{
// Assignment
   if (&other != this) {
      fEvent = other.fEvent;
      fEvslot = other.fEvslot;
      fParticle = other.fParticle;
      fPDG = other.fPDG;
      fG5code = other.fG5code;
      fCharge = other.fCharge;
      fProcess = other.fProcess;
      fIzero = other.fIzero;
      fNsteps = other.fNsteps;
      fSpecies = other.fSpecies;
      fStatus = other.fStatus;
      fMass = other.fMass;
      fXpos = other.fXpos;
      fYpos = other.fYpos;
      fZpos = other.fZpos;
      fXdir = other.fXdir;
      fYdir = other.fYdir;
      fZdir = other.fZdir;
      fP = other.fP;
      fE = other.fE;
      fPstep = other.fPstep;
      fStep = other.fStep;
      fSnext = other.fSnext;
      fSafety = other.fSafety; 
      fFrombdr = other.fFrombdr;
      fPending = other.fPending;
      fPath = new TGeoBranchArray(*other.fPath);
      fNextpath = new TGeoBranchArray(*other.fNextpath);
   }
   return *this;
}
   
//______________________________________________________________________________
GeantTrack::~GeantTrack()
{
// Destructor.
   delete path;
   delete nextpath;
}   

//______________________________________________________________________________
void GeantTrack::ReadFromVector(const GeantTrack_v arr. Int_t i)
{
// Fill track from array
   fEvent = arr.fEventV[i];
   fEvslot = arr.fEvslotV[i];
   fParticle = arr.fParticleV[i];
   fPDG = arr.fPDGV[i];
   fG5code = arr.fG5codeV[i];
   fCharge = arr.fChargeV[i];
   fProcess = arr.fProcessV[i];
   fIzero = arr.fIzeroV[i];
   fNsteps = arr.fNstepsV[i];
   fSpecies = arr.fSpeciesV[i];
   fStatus = arr.fStatusV[i];
   fMass = arr.fMassV[i];
   fXpos = arr.fXposV[i];
   fYpos= arr.fYposV[i];
   fZpos = arr.fZposV[i];
   fXdir = arr.fXdirV[i];
   fYdir = arr.fYdirV[i];
   fZdir = arr.fZdirV[i];
   fP = arr.fPV[i];
   fE = arr.fEV[i];
   fPstep = arr.fPstepV[i];
   fStep = arr.fStepV[i];
   fSnext = arr.fSnextV[i];
   fSafety = arr.fSafetyV[i];
   fFrombdr = arr.fFrombdrV[i];
   fPending = arr.fPendingV[i];
   fPath = arr.fPathV[i];
   fNextpath = arr.fNextpathV[i];
}

//______________________________________________________________________________
void GeantTrack::Reset()
{
// Resets track content.
   fEvent = -1;
   fEvslot = -1;
   fParticle = -1;
   fPDG = 0;
   fG5code = 0;
   fCharge = 0;
   fProcess = -1;
   fIzero = 0;
   fNsteps = 0;
   fSpecies = kHadron;
   fStatus = kAlive;
   fMass = 0.;
   fXpos = 0.;
   fYpos = 0.;
   fZpos = 0.;
   fXdir = 0.;
   fYdir = 0.;
   fZdir = 0.;
   fP = 0.;
   fE = 0.;
   fPstep = 1.E20;
   fStep = 0.;
   fSnext = 0.;
   fSafety = 0.;
   fFrombdr = false;
   fPending = false;
}   

//______________________________________________________________________________
void GeantTrack::Print(Int_t) const {
   TString spath;
//   if (path) path->GetPath(spath);
   Printf("=== Track %d (ev=%d): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Dir:(%f,%f,%f) P:%g E:%g snext=%g safety=%g nsteps=%d",
           fParticle,fEvent, fProcess,fPstep,fCharge,fXpos,fYpos,fZpos,fXdir,fYdir,fZdir,P(),fE,fSnext,fSafety,fNsteps);
}


ClassImp(GeantTrack_v)

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v()
             :TObject(),
              fNtracks(0),fMaxtracks(0),fNselected(0),fHoles(),fSelected(),fCompact(true),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
// Dummy ctor.
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(Int_t size)
             :TObject(),
              fNtracks(0),fMaxtracks(0),fNselected(0),fHoles(size),fSelected(size),fCompact(true),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
// Constructor with maximum capacity.
   Resize(size);
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(const GeantTrack_v &track_v)
             :TObject(),
              fNtracks(track_v.fNtracks),fMaxtracks(track_v.fMaxTracks),fNselected(track_v.fNselected),fHoles(track_v.fHoles),
              fSelected(track_v.fSelected),fCompact(track_v.fCompact),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
// Copy constructor
   Int_t size = track_v.fMaxtracks;
   fBuf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
   memcpy(fBuf, track_v.fBuf, fNtracks*sizeof(GeantTrack));
   AssignInBuffer(fBuf, size);
}   

//______________________________________________________________________________
GeantTrack_v &GeantTrack_v()::operator=(const GeantTrack_v &track_v)
{
// Assignment operator
   if (&other != this) {
      Int_t size = track_v.fMaxtracks;
      if (fMaxtracks<size) {
         delete [] fBuf;
         fBuf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
      }
      fNtracks = track_v.fNtracks;
      fMaxtracks = size;
      fNselected = track_v.fNselected;
      fHoles = track_v.fHoles;
      fSelected = track_v.fSelected;
      fCompact = track_v.fCompact();
      memcpy(fBuf, track_v.fBuf, size*sizeof(GeantTrack));
      AssignInBuffer(fBuf, size);
   }
   return *this;   
}   

//______________________________________________________________________________
GeantTrack_v()::~GeantTrack_v()
{
// Destructor.
   _mm_free(fBuf);
}
   
//______________________________________________________________________________  
void GeantTrack_v::AssignInBuffer(const char *buff, Int_t size)
{
// Assign all internal class arrays in the supplied buffer, padded by supplied 
// size.
   const Int_t size_intn = size*sizeof(Int_t);
   const Int_t size_doublen = size*sizeof(Double_t);
   const Int_t size_booln = size*sizeof(Bool_t);
   char *buf = (char*)buff;
   fEventV = (Int_t*)buf;
   buf += size_intn;
   fEvslotV = (Int_t*)buf;
   buf += size_intn;
   fParticleV = (Int_t*)buf;
   buf += size_intn;
   fPDGV = (Int_t*)buf;
   buf += size_intn;
   fG5codeV = (Int_t*)buf;
   buf += size_intn;
   fChargeV = (Int_t*)buf;
   buf += size_intn;
   fProcessV = (Int_t*)buf;
   buf += size_intn;
   fIzeroV = (Int_t*)buf;
   buf += size_intn;
   fNstepsV = (Int_t*)buf;
   buf += size_intn;
   fSpeciesV = (Species_t*)buf;
   buf += size*sizeof(Species_t);
   fStatusV = (TrackStatus_t*)buf;
   buf += size*sizeof(TrackStatus_t);
   fMassV = (Double_t*)buf;
   buf += size_doublen;
   fXposV = (Double_t*)buf;
   buf += size_doublen;
   fYposV = (Double_t*)buf;
   buf += size_doublen;
   fZposV = (Double_t*)buf;
   buf += size_doublen;
   fXdirV = (Double_t*)buf;
   buf += size_doublen;
   fYdirV = (Double_t*)buf;
   buf += size_doublen;
   fZdirV = (Double_t*)buf;
   buf += size_doublen;
   fPV = (Double_t*)buf;
   buf += size_doublen;
   fEV = (Double_t*)buf;
   buf += size_doublen;
   fPstepV = (Double_t*)buf;
   buf += size_doublen;
   fStepV = (Double_t*)buf;
   buf += size_doublen;
   fSnextV = (Double_t*)buf;
   buf += size_doublen;
   fSafetyV = (Double_t*)buf;
   buf += size_doublen;
   fFrombdrV = (Bool_t*)buf;
   buf += size_booln;
   fPendingV = (Bool_t*)buf;
   buf += size_booln;
   fPathV = (TGeoBranchArray**)buf;
   buf += size*sizeof(TGeoBranchArray*);
   fNextpathV = (TGeoBranchArray**)buf;
   buf += size*sizeof(TGeoBranchArray*);
}

//______________________________________________________________________________  
void GeantTrack_v::CopyToBuffer(const char *buff, Int_t size)
{
// Copy existing track arrays into new buffer, padded by supplied size
   const Int_t size_int = fNtracks*sizeof(Int_t);
   const Int_t size_double = fNtracks*sizeof(Double_t);
   const Int_t size_intn = size*sizeof(Int_t);
   const Int_t size_doublen = size*sizeof(Double_t);
   char *buf = (char*)buff;
   fEventV = (Int_t*)buf;   
   memcpy(buf, fEventV, size_int);
   buf += size_intn;
   fEventV = (Int_t*)buf;
   memcpy(buf, fEvslotV, size_int);
   buf += size_intn;
   fEvslotV = (Int_t*)buf;
   memcpy(buf, fParticleV, size_int);
   buf += size_intn;
   fParticleV = (Int_t*)buf;
   memcpy(buf, fPDGV, size_int);
   buf += size_intn;
   fPDGV = (Int_t*)buf;
   memcpy(buf, fG5codeV, size_int);
   buf += size_intn;
   fG5codeV = (Int_t*)buf;
   memcpy(buf, fChargeV, size_int);
   buf += size_intn;
   fChargeV = (Int_t*)buf;
   memcpy(buf, fProcessV, size_int);
   buf += size_intn;
   fProcessV = (Int_t*)buf;
   memcpy(buf, fIzeroV, size_int);
   buf += size_intn;
   fIzeroV = (Int_t*)buf;
   memcpy(buf, fNstepsV, size_int);
   buf += size_intn;
   fNstepsV = (Int_t*)buf;
   memcpy(buf, fSpeciesV, fNtracks*sizeof(Species_t));
   buf += size*sizeof(Species_t);
   fSpeciesV = (Species_t*)buf;
   memcpy(buf, fStatusV, fNtracks*sizeof(TrackStatus_t));
   buf += size*sizeof(TrackStatus_t);
   fStatusV = (TrackStatus_t*)buf;
   memcpy(buf, fMassV, size_double);
   buf += size_doublen;
   fMassV = (Double_t*)buf;
   memcpy(buf, fXposV, size_double);
   buf += size_doublen;
   fXposV = (Double_t*)buf;
   memcpy(buf, fYposV, size_double);
   buf += size_doublen;
   fYposV = (Double_t*)buf;
   memcpy(buf, fZposV, size_double);
   buf += size_doublen;
   fZposV = (Double_t*)buf;
   memcpy(buf, fXdirV, size_double);
   buf += size_doublen;
   fXdirV = (Double_t*)buf;
   memcpy(buf, fYdirV, size_double);
   buf += size_doublen;
   fYdirV = (Double_t*)buf;
   memcpy(buf, fZdirV, size_double);
   buf += size_doublen;
   fZdirV = (Double_t*)buf;
   memcpy(buf, fPV, size_double);
   buf += size_doublen;
   fPV = (Double_t*)buf;
   memcpy(buf, fEV, size_double);
   buf += size_doublen;
   fEV = (Double_t*)buf;
   memcpy(buf, fPstepV, size_double);
   buf += size_doublen;
   fPstepV = (Double_t*)buf;
   memcpy(buf, fStepV, size_double);
   buf += size_doublen;
   fStepV = (Double_t*)buf;
   memcpy(buf, fSnextV, size_double);
   buf += size_doublen;
   fSnextV = (Double_t*)buf;
   memcpy(buf, fSafetyV, size_double);
   buf += size_doublen;
   fSafetyV = (Double_t*)buf;
   memcpy(buf, fFrombdrV, fNtracks*sizeof(Bool_t));
   buf += size*sizeof(Bool_t);
   fFrombdrV = (Bool_t*)buff;
   memcpy(buf, fPendingV, fNtracks*sizeof(Bool_t));
   buf += size*sizeof(Bool_t);
   fPendingV = (Bool_t*)buff;
   memcpy(buf, fPathV, fNtracks*sizeof(TGeoBranchArray*));
   buf += size*sizeof(TGeoBranchArray*);
   fPathV = (TGeoBranchArray**)buf;
   memcpy(buf, fNextpathV, fNtracks*sizeof(TGeoBranchArray*));
   buf += size*sizeof(TGeoBranchArray*);
   fNextpathV = (TGeoBranchArray**)buf;
}

//______________________________________________________________________________  
void GeantTrack_v::Resize(Int_t newsize)
{
// Resize the container.
   Int_t size = round_up_align(newsize);
   if (size<fNtracks) {
      printf("Error: Cannot resize to less than current track content\n");
      return;
   } 
   fHoles.SetBitNumber(size-1, kFALSE);        
   char* buf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
   if (!fBuf) {
      // All arrays are contiguous in a single buffer and aligned with the
      // same padding ALIGN_PADDING
      fBuf = buf;
      fMaxtracks = size;
      AssignInBuffer(buf, size);
   } else {
      // Resize container
      CopyToBuffer(buf, size);
      delete [] fBuf;
      fBuf = buf;
      fMaxtracks = size;
   }   
}

//______________________________________________________________________________  
void GeantTrack_v::AddTrack(const GeantTrack &track)
{
// Add new track to the array. If addition is done on top of non-compact array, 
// the track will be inserted without updating the number of tracks. 
   Int_t itrack = fNtracks;
   if (!fCompact) itrack = fHoles.FirstSetBit();
   if (itrack==fMaxtracks) Resize(2*fMaxtracks);
   fEventV[itrack] = track.Event();
   fEvslotV[itrack] = track.EventSlot();
   fParticleV[itrack] = track.Particle();
   fPDGV[itrack] = track.PDG();
   fG5codeV[itrack] = track.G5code();
   fProcessV[itrack] = track.Process();
   fIzeroV[itrack] = track.Izero();
   fNstepsV[itrack] = track.GetNsteps();
   fSpeciesV[itrack] = track.Species();
   fStatusV[itrack] = track.Status();
   fMassV[itrack] = track.Mass();
   fXposV[itrack] = track.PosX();
   fYposV[itrack] = track.PosY();
   fZposV[itrack] = track.PosZ();
   fXdirV[itrack] = track.DirX();
   fYdirV[itrack] = track.DirY();
   fZdirV[itrack] = track.DirZ();
   fPV[itrack] = track.P();
   fEV[itrack] = track.E();
   fPstepV[itrack] = track.GetPstep();
   fStepV[itrack] = track.GetStep();
   fSnextV[itrack] = track.GetSnext();
   fSafetyV[itrack] = track.GetSafety();
   fFrombdrV[itrack] = track.FromBoundary();
   fPendingV[itrack] = track.Pending();
   fPathV[itrack] = track.GetPath();
   fNextpathV[itrack] = track.GetNextPath();
   if (itrack==fNtracks) fNtracks++;
}   

//______________________________________________________________________________  
void GeantTrack_v::GetTrack(Int_t i, GeantTrack &track) {
// Extract a single track from array.
   track.ReadFromVector(this);
}
   
//______________________________________________________________________________  
void GeantTrack_v::AddTrack(const GeantTrack_v &arr, Int_t i)
{
// Add track from different array
   if (fNtracks==fMaxtracks) Resize(2*fMaxtracks);
   
   fEventV[fNtracks] = arr.fEventV[i];
   fEvslotV[fNtracks] = arr.fEvslotV[i];
   fParticleV[fNtracks] = arr.fParticleV[i];
   fPDGV[fNtracks] = arr.fPDGV[i];
   fG5codeV[fNtracks] = arr.fG5codeV[i];
   fChargeV[fNtracks] = arr.fChargeV[i];
   fProcessV[fNtracks] = arr.fProcessV[i];
   fIzeroV[fNtracks] = arr.fIzeroV[i];
   fNstepsV[fNtracks] = arr.fNstepsV[i];
   fSpeciesV[fNtracks] = arr.fSpeciesV[i];
   fStatusV[fNtracks] = arr.fStatusV[i];
   fMassV[fNtracks] = arr.fMassV[i];
   fXposV[fNtracks] = arr.fXposV[i];
   fYposV[fNtracks] = arr.fYposV[i];
   fZposV[fNtracks] = arr.fZposV[i];
   fXdirV[fNtracks] = arr.fXdirV[i];
   fYdirV[fNtracks] = arr.fYdirV[i];
   fZdirV[fNtracks] = arr.fZdirV[i];
   fPV[fNtracks] = arr.fPV[i];
   fEV[fNtracks] = arr.fEV[i];
   fPstepV[fNtracks] = arr.fPstepV[i];
   fStepV[fNtracks] = arr.fStepV[i];
   fSnextV[fNtracks] = arr.fSnextV[i];
   fSafetyV[fNtracks] = arr.fSafetyV[i];
   fFrombdrV[fNtracks] = arr.fFrombdrV[i];
   fPendingV[fNtracks] = arr.fPendingV[i];
   fPathV[fNtracks] = arr.fPathV[i];
   fNextpathV[fNtracks] = arr.fNextpathV[i];

   fSelected.ResetBitNumber(fNtracks);
   fHoles.ResetBitNumber(fNtracks);
   fNtracks++;
}

//______________________________________________________________________________  
void GeantTrack_v::AddTracks(const GeantTrack_v &arr, Int_t istart, Int_t iend)
{
// Add tracks from different array
   Int_t ncpy = iend-istart+1;
   if (fNtracks+ncpy>=fMaxtracks) Resize(TMath::Max(2*fMaxtracks, fNtracks+ncpy));   
   memcpy(&fEventV[fNtracks], &arr.fEventV[istart], ncpy*sizeof(Int_t));
   memcpy(&fEvslotV[fNtracks], &arr.fEvslotV[i], ncpy*sizeof(Int_t));
   memcpy(&fParticleV[fNtracks], &arr.fParticleV[i], ncpy*sizeof(Int_t));
   memcpy(&fPDGV[fNtracks], &arr.fPDGV[i], ncpy*sizeof(Int_t));
   memcpy(&fG5codeV[fNtracks], &arr.fG5codeV[i], ncpy*sizeof(Int_t));
   memcpy(&fChargeV[fNtracks], &arr.fChargeV[i], ncpy*sizeof(Int_t));
   memcpy(&fProcessV[fNtracks], &arr.fProcessV[i], ncpy*sizeof(Int_t));
   memcpy(&fIzeroV[fNtracks], &arr.fIzeroV[i], ncpy*sizeof(Int_t));
   memcpy(&fNstepsV[fNtracks], &arr.fNstepsV[i], ncpy*sizeof(Int_t));
   memcpy(&fSpeciesV[fNtracks], &arr.fSpeciesV[i], ncpy*sizeof(Species_t));
   memcpy(&fStatusV[fNtracks], &arr.fStatusV[i], ncpy*sizeof(Trackstatus_t));
   memcpy(&fMassV[fNtracks], &arr.fMassV[i], ncpy*sizeof(Double_t));
   memcpy(&fXposV[fNtracks], &arr.fXposV[i], ncpy*sizeof(Double_t));
   memcpy(&fYposV[fNtracks], &arr.fYposV[i], ncpy*sizeof(Double_t));
   memcpy(&fZposV[fNtracks], &arr.fZposV[i], ncpy*sizeof(Double_t));
   memcpy(&fXdirV[fNtracks], &arr.fXdirV[i], ncpy*sizeof(Double_t));
   memcpy(&fYdirV[fNtracks], &arr.fYdirV[i], ncpy*sizeof(Double_t));
   memcpy(&fZdirV[fNtracks], &arr.fZdirV[i], ncpy*sizeof(Double_t));
   memcpy(&fPV[fNtracks], &arr.fPV[i], ncpy*sizeof(Double_t));
   memcpy(&fEV[fNtracks], &arr.fEV[i], ncpy*sizeof(Double_t));
   memcpy(&fPstepV[fNtracks], &arr.fPstepV[i], ncpy*sizeof(Double_t));
   memcpy(&fStepV[fNtracks], &arr.fStepV[i], ncpy*sizeof(Double_t));
   memcpy(&fSnextV[fNtracks], &arr.fSnextV[i], ncpy*sizeof(Double_t));
   memcpy(&fSafetyV[fNtracks], &arr.fSafetyV[i], ncpy*sizeof(Double_t));
   memcpy(&fFrombdrV[fNtracks], &arr.fFrombdrV[i], ncpy*sizeof(Bool_t));
   memcpy(&fPendingV[fNtracks], &arr.fPendingV[i], ncpy*sizeof(Bool_t));
   memcpy(&fPathV[fNtracks], &arr.fPathV[i], ncpy*sizeof(TGeoBranchArray*));
   memcpy(&fNextpathV[fNtracks], &arr.fNextpathV[i], ncpy*sizeof(TGeoBranchArray*));
   fSelected.ResetBitNumber(fNtracks+ncpy-1);
   fHoles.ResetBitNumber(fNtracks+ncpy-1);
   fNtracks += ncpy;
}

//______________________________________________________________________________  
void GeantTrack_v::SwapTracks(Int_t i, Int_t j)
{
// Swap two tracks in the container
   Double_t tdbl;
   Int_t    tint;
   Bool_t   tbool;
   TGeoBranchArray *tptr;
   tint = fEventV[i]; fEventV[i] = fEventV[j]; fEventV[j] = tint;
   tint = fEvslotV[i]; fEvslotV[i] = fEvslotV[j]; fEvslotV[j] = tint;
   tint = fParticleV[i]; fParticleV[i] = fParticleV[j]; fParticleV[j] = tint;
   tint = fPDGV[i]; fPDGV[i] = fPDGV[j]; fPDGV[j] = tint;
   tint = fG5codeV[i]; fG5codeV[i] = fG5codeV[j]; fG5codeV[j] = tint;
   tint = fChargeV[i]; fChargeV[i] = fChargeV[j]; fChargeV[j] = tint;
   tint = fProcessV[i]; fProcessV[i] = fProcessV[j]; fProcessV[j] = tint;
   tint = fIzeroV[i]; fIzeroV[i] = fIzeroV[j]; fIzeroV[j] = tint;
   tint = fNstepsV[i]; fNstepsV[i] = fNstepsV[j]; fNstepsV[j] = tint;
   tint = fSpeciesV[i]; fSpeciesV[i] = fSpeciesV[j]; fSpeciesV[j] = tint;
   tint = fStatusV[i]; fStatusV[i] = fStatusV[j]; fStatusV[j] = tint;
   tdbl = fMassV[i]; fMassV[i] = fMassV[j]; fMassV[j] = tdbl;
   tdbl = fMassV[i]; fMassV[i] = fMassV[j]; fMassV[j] = tdbl;
   tdbl = fXposV[i]; fXposV[i] = fXposV[j]; fXposV[j] = tdbl;
   tdbl = fYposV[i]; fYposV[i] = fYposV[j]; fYposV[j] = tdbl;
   tdbl = fZposV[i]; fZposV[i] = fZposV[j]; fZposV[j] = tdbl;
   tdbl = fXdirV[i]; fXdirV[i] = fXdirV[j]; fXdirV[j] = tdbl;rent vector-based steering framework for particle 
   tdbl = fYdirV[i]; fYdirV[i] = fYdirV[j]; fYdirV[j] = tdbl;
   tdbl = fZdirV[i]; fZdirV[i] = fZdirV[j]; fZdirV[j] = tdbl;
   tdbl = fPV[i]; fPV[i] = fPV[j]; fPV[j] = tdbl;
   tdbl = fEV[i]; fEV[i] = fEV[j]; fEV[j] = tdbl;
   tdbl = fPstepV[i]; fPstepV[i] = fPstepV[j]; fPstepV[j] = tdbl;
   tdbl = fSnextV[i]; fSnextV[i] = fSnextV[j]; fSnextV[j] = tdbl;
   tdbl = fSafetyV[i]; fSafetyV[i] = fSafetyV[j]; fSafetyV[j] = tdbl;
   tbool = fFrombdrV[i]; fFrombdrV[i] = fFrombdrV[j]; fFrombdrV[j] = tbool;
   tbool = fPendingV[i]; fPendingV[i] = fPendingV[j]; fPendingV[j] = tbool;
   tptr = fPathV[i]; fPathV[i] = fPathV[j]; fPathV[j] = tptr;
   tptr = fNextpathV[i]; fNextpathV[i] = fNextpathV[j]; fNextpathV[j] = tptr;
   Bool_t sel = fSelected.TestBitNumber(j);
   fSelected.SetBitNumber(j, fSelected.TestBitNumber(i));
   fSelected.SetBitNumber(i, sel);
}

//______________________________________________________________________________  
void GeantTrack_v::ReplaceTrack(Int_t i, Int_t j)
// Replace content of track i with the one of track j
   fEventV[i] = fEventV[j];
   fEvslotV[i] = fEvslotV[j];
   fParticleV[i] = fParticleV[j];
   fPDGV[i] = fPDGV[j];
   fG5codeV[i] = fG5codeV[j];
   fChargeV[i] = fChargeV[j];
   fProcessV[i] = fProcessV[j];
   fIzeroV[i] = fIzeroV[j];
   fNstepsV[i] = fNstepsV[j];
   fSpeciesV[i] = fSpeciesV[j];
   fStatusV[i] = fStatusV[j];
   fMassV[i] = fMassV[j];
   fMassV[i] = fMassV[j];
   fXposV[i] = fXposV[j];
   fYposV[i] = fYposV[j];
   fZposV[i] = fZposV[j];
   fXdirV[i] = fXdirV[j];
   fYdirV[i] = fYdirV[j];
   fZdirV[i] = fZdirV[j];
   fPV[i] = fPV[j];
   fEV[i] = fEV[j];
   fPstepV[i] = fPstepV[j];
   fSnextV[i] = fSnextV[j];
   fSafetyV[i] = fSafetyV[j];
   fFrombdrV[i] = fFrombdrV[j];
   fPendingV[i] = fPendingV[j];
   fPathV[i] = fPathV[j];
   fNextpathV[i] = fNextpathV[j];
   fSelected.SetBitNumber(i, fSelected.TestBitNumber(j));
}   

//______________________________________________________________________________  
void GeantTrack_v::RemoveTracks(Int_t from, Int_t to)
{
// Remove tracks from the container
   Int_t ncpy = fNtracks-to;
   memmove(&fEventV[from], &fEventV[to+1], ncpy*sizeof(Int_t));
   memmove(&fEvslotV[from], &fEvslotV[to+1], ncpy*sizeof(Int_t));
   memmove(&fParticleV[from], &fParticleV[to+1], ncpy*sizeof(Int_t));
   memmove(&fPDGV[from], &fPDGV[to+1], ncpy*sizeof(Int_t));
   memmove(&fG5codeV[from], &fG5codeV[to+1], ncpy*sizeof(Int_t));
   memmove(&fChargeV[from], &fChargeV[to+1], ncpy*sizeof(Int_t));
   memmove(&fProcessV[from], &fProcessV[to+1], ncpy*sizeof(Int_t));
   memmove(&fIzeroV[from], &fIzeroV[to+1], ncpy*sizeof(Int_t));
   memmove(&fNstepsV[from], &fNstepsV[to+1], ncpy*sizeof(Int_t));
   memmove(&fSpeciesV[from], &fSpeciesV[to+1], ncpy*sizeof(Species_t));
   memmove(&fStatusV[from], &fStatusV[to+1], ncpy*sizeof(Trackstatus_t));
   memmove(&fMassV[from], &fMassV[to+1], ncpy*sizeof(Double_t));
   memmove(&fXposV[from], &fXposV[to+1], ncpy*sizeof(Double_t));
   memmove(&fYposV[from], &fYposV[to+1], ncpy*sizeof(Double_t));
   memmove(&fZposV[from], &fZposV[to+1], ncpy*sizeof(Double_t));
   memmove(&fXdirV[from], &fXdirV[to+1], ncpy*sizeof(Double_t));
   memmove(&fYdirV[from], &fYdirV[to+1], ncpy*sizeof(Double_t));
   memmove(&fZdirV[from], &fZdirV[to+1], ncpy*sizeof(Double_t));
   memmove(&fPV[from], &fPV[to+1], ncpy*sizeof(Double_t));
   memmove(&fEV[from], &fEV[to+1], ncpy*sizeof(Double_t));
   memmove(&fPstepV[from], &fPstepV[to+1], ncpy*sizeof(Double_t));
   memmove(&fStepV[from], &fStepV[to+1], ncpy*sizeof(Double_t));
   memmove(&fSnextV[from], &fSnextV[to+1], ncpy*sizeof(Double_t));
   memmove(&fSafetyV[from], &fSafetyV[to+1], ncpy*sizeof(Double_t));
   memmove(&fFrombdrV[from], &fFrombdrV[to+1], ncpy*sizeof(Bool_t));
   memmove(&fPendingV[from], &fPendingV[to+1], ncpy*sizeof(Bool_t));
   memmove(&fPathV[from], &fPathV[to+1], ncpy*sizeof(TGeoBranchArray*));
   memmove(&fNextpathV[from], &fNextpathV[to+1], ncpy*sizeof(TGeoBranchArray*));
   fNtracks -= ncpy;
   fSelected.ResetAllBits();
   fNselected = 0;
}   
   
//______________________________________________________________________________  
Int_t GeantTrack_v::Compact(GeantTrack_v *moveto)
{
// Compact the holes in the array. Return number of active elements. This will
// lose the track fields in the holes, so information from the holes has to be
// copied beforehand
   fCompact = kTRUE;
   if (fNtracks == 0) return 0;
   Int_t firsthole = fHoles.FirstSetBit();
   while (firsthole<fNtracks) {
      Int_t lastactive = fHoles.LastNullBit(fNtracks-1);
      if (lastactive < fNtracks) {
         fNtracks = lastactive+1;
         if (firsthole==fNtracks) return fNtracks;
      } else {
         // No active tracks left
         fNtracks = 0;
         return 0;
      }
      // replace content of first hole with the last active track
      if (moveto) moveto->AddTrack(*this, firsthole);
      ReplaceTrack(firsthole, lastactive);
      fHoles.SetBitNumber(firsthole, false);
      fHoles.SetBitNumber(lastactive, true);
      firsthole = fHoles.FirstSetBit(firsthole+1);
      fNtracks--;
   }
   fSelected->ResetAllBits();
   fNselected = 0;
   return fNtracks;
}

//______________________________________________________________________________  
Int_t GeantTrack_v::Reshuffle()
{
// Reshuffle tracks according the selection mask. The selected tracks will be
// moved at the beginning of the array. Tracks should be compacted before.
   if (fNtracks == 0) return 0;
   fNselected = fNtracks;
   Int_t firsthole = fSelected.FirstNullBit();
   while (firsthole<fNselected) {
      Int_t lastsel = fSelected.LastSetBit(fNselected-1);
      if (lastsel >= fNselected) return 0;
      fNselected = lastsel+1;
      if (firsthole==fNselected) return fNselected;
      // exchange tracks pointed by firsthole and lastactive
      SwapTracks(firsthole, lastsel);
      fSelected.SetBitNumber(firsthole, true);
      fSelected.SetBitNumber(lastactive, false);
      firsthole = fSelected.FirstNullBit(firsthole+1);
      fNselected--;
   }
   return fNselected;
}
      
//______________________________________________________________________________
Int_t GeantTrack_v::PropagateStraight(Int_t ntracks, Double_t *crtstep)
{
// Propagate first ntracks along a straight line (neutral particles, no mag. 
// field or for last tiny step). The array has to be reshuffled after selecting
// the interesting particles using Select method.
// The crossing tracks get masked as holes in the array.

// Find next volume
   Int_t icrossed = 0;
   for (Int_t i=0; i<ntracks; i++) {
      if (fFrombdrV[i]) {
         *fPathV[i] = *fNextPathV[i];
         fStatus[i] = kBoundary;
//         MarkRemoved(i);
         icrossed++;
      } else {
         
   }   
   for (Int_t i=0; i<ntracks; i++) {
      fPstep[i] -= crtstep[i];
      fSafety[i] = 0;
      // Change path to reflect the physical volume for the current track; The
      // relevant path is fPath[i] if the frombdr flag is not set or fNextpath[i]
      // otherwise
      fXpos[i] += crtstep[i]*fXdir[i];
      fYpos[i] += crtstep[i]*fYdir[i];
      fZpos[i] += crtstep[i]*fZdir[i];
      fNsteps[i]++;
   }   
   return icrossed;
}   

//______________________________________________________________________________
void GeantTrack_v::PropagateInVolume(const Double_t *crtstep)
{
// Propagate the selected tracks with crtstep values. The method is to be called
// only with  charged tracks in magnetic field. The tracks array has to be reshuffled.
   Double_t c = 0.;
   Double_t ptot = 0;
   const Double_t *point = 0;
   const Double_t *newdir = 0;
   const Double_t bmag = gPropagator->fBmag;
   GeantThreadData *td = gPropagator->fThreadData[TGeoManager::ThreadId();
   TGeoHelix *fieldp = td->fFieldPropagator;
   for (Int_t i=0; i<fNtracks; i++) {
      // Reset relevant variables
      fFrombdr[i] = kFALSE;
      fPstep[i] -= crtstep[i];
      fSafety[i] -= crtstep[i];
      if (fSafety[i]<0.) fSafety = 0.;
      fStep[i] += crtstep[i];
      // Set curvature, charge
      c = std::fabs(kB2C*bmag/Pt(i));
// NOTE: vectorized treatment in TGeoHelix
      fieldp->SetXYcurvature(c);
      fieldp->SetCharge(fCharge);
      fieldp->SetHelixStep(std::fabs(TMath::TwoPi()*Pz(i)/(c*Pt(i))));
      fieldp->InitPoint(fXpos[i],fYpos[i],fZpos[i]);
      fieldp->InitDirection(fXdir[i],fYdir[i], fZdir[i]);
      fieldp->UpdateHelix();
      fieldp->Step(crtstep[i]);
      point = fieldp->GetCurrentPoint();
      newdir = fieldp->GetCurrentDirection();
      fXpos[i] = point[0]; fYpos[i] = point[1]; fZpos[i] = point[2];
      fXdir[i] = newdir[0]; fYdir[i] = newdir[1]; fZdir[i] = newdir[2];
   }
}   

//______________________________________________________________________________
void GeantTrack_v::NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep, 
                       const Double_t *x, const Double_t *y, const Double_t *z,
                       const Double_t *dirx, const Double_t *diry, const Double_t *dirz,
                       TGeoBranchArray **pathin, TGeoBranchArray **pathout, 
                       Double_t *step, Double_t *safe, Bool_t *isonbdr)
{
// Vector version of TGeo FNB (To be implemented the vectorized navigator)
// Apply to all particles (charged or not).
//    pstep = proposed steps by physics
//    x,y,z, dirx,diry, dirz = initial positions and directions
//    safety = safety values for the initial points
//    step = distances to next boundary (or to physics step if closer)
//    isonbdr = starting points on boundary flags (used to decide on safety computation)
//              use also as output to notify if the step is boundary or physics
//    pathin = starting paths
//    pathout = final path after propagation to next boundary
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();   
   for (Int_t i=0; i<ntracks; i++) {
      nav->ResetState();
      nav->SetCurrentPoint(x[i], y[i], z[i]);
      nav->SetCurrentDirection(dirx[i], diry[i], dirz[i]);
      pathin[i]->UpdateNavigator(nav);
      nav->SetLastSafetyForPoint(safe[i], x[i], y[i], z[i]);
      nav->FindNextBoundaryAndStep(TMath::Min(1.E20, pstep[i]), isonbdr[i]);
      step[i] = TMath::Max(2*gTolerance,nav->GetStep());
      safe[i] = nav->GetSafeDistance();
      pathout[i]->InitFromNavigator(nav);
      isonbdr[i] = nav->IsOnBoundary();
   }      
}

//______________________________________________________________________________
void GeantTrack_v::NavIsSameLocation(TGeoBranchArray *start, TGeoBranchArray *end, Bool_t *same)
{
// Implementation of TGeoNavigator::IsSameLocation with vector input
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   for (Int_t i=0; i<fNtracks; i++) {   
      nav->ResetState();
      nav->SetCurrentPoint(fXpos[i], fYpos[i], fZpos[i]);
      nav->SetCurrentDirection(fXdir[i], fYdir[i], fZdir[i]);
      start[i]->UpdateNavigator(nav);
      same[i] = nav->IsSameLocation(fXpos[i],fYpos[i],fZpos[i],kTRUE);
      if (same[i]) end->InitFromNavigator(nav);
   }
}   
//______________________________________________________________________________
void GeantTrack_v::PropagateBack(Int_t itr)
{
// This method is to be called after a successful crossing made by PropagateInField
// to put the particle at the entry point.
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   fNextPathV[i]->UpdateNavigator(nav);
   TGeoNode *checked = nav->GetCurrentNode();
   TGeoVolume *vol = checked->GetVolume();
   Double_t ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;
   Bool_t outside = nav->IsOutside();
   // Swap track direction and compute distance back to boundary
   fXdir[i] = -newdir[0]; fYdir[i] = -newdir[1]; fZdir[i] = -newdir[2];
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   TGeoNode *node1 = 0;
   TGeoNode *node2 = 0;
   if (level < fPathV[i]->GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         node1 = nav->GetMother(level-lev);
         node2 = fPath[i]->GetNode(lev);
         if (node1 == node2) {
            if (lev==level) entering = kFALSE;
         } else {
         // different nodes at some level -> entering current node
            break;
         }   
      }
   }   
   if (!entering) {
      checked = fPath[i]->GetNode(level+1);
      if (!checked) return 0;
   }   
   nav->MasterToLocal(&fXpos[i], local);
   nav->MasterToLocalVect(&fXdir[i], ldir);
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }   
   if (useDebug && (debugTrk<0 || itr==debugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }   
   if (delta>crtstep) {
      if (useDebug && (debugTrk<0 || itr==debugTrk)) {
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
   fXpos[i] += delta*fXdir[i];
   fYpos[i] += delta*fYdir[i];
   fZpos[i] += delta*fZdir[i];
   
//______________________________________________________________________________
Int_t GeantTrack_v::PropagateInField(Int_t ntracks, const Double_t *crtstep)
{
// Propagate with crtstep using the helix propagator. Mark crossing tracks as holes.
   Int_t icrossed = 0;
//   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
//   Int_t tid = nav->GetThreadId();
//   GeantThreadData *td = gPropagator->fThreadData[tid];
//   Bool_t useDebug = gPropagator->fUseDebug;
//   Int_t debugTrk = gPropagator->fDebugTrk;
   Bool_t *same = new Bool_t[ntracks];
   PropagateInVolume(crtstep);
   NavIsSameLocation(fPathV, fNextPathV, same);
   for (Int_t itr=0; itr<ntracks; itr++) {
      if (same[itr]) continue;      
      // Boundary crossed
      PropagateBack(itr);
   // Create a new branch array
   fPath[i]->InitFromNavigator(nav);
//   if (vol->IsAssembly()) Printf("### ERROR ### Entered assembly %s", vol->GetName());
   MarkRemoved(i);
   icrossed++;
   return icrossed;
}   

//______________________________________________________________________________
void GeantTrack_v::ComputeTransportLength(Int_t ntracks)
{
// Computes snext and safety for an array of tracks. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.
   static Int_t icalls = 0;
   icalls++;  
   Int_t itr;
   Bool_t isOnBoundary = kFALSE;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   NavFindNextBoundaryAndStep(ntracks, fPstepV, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV,
                              fPathV, fNextPathV, fStepV, fSafeV, fFrombdrV);
   for (itr=0; itr<ntracks; itr++) {
      if (fNextPathV[i]->IsOutside() || fSnextV[i]>1.E19) fStatusV[i] = kExitingSetup;
      if (fFrombdrV[i] && fSnextV[i]<2.*gTolerance) {
         // Make sure track crossed
         fIzeroV[i]++;
         if (fIzeroV[i] > 10) {
            fStatus[i] = kKilled;
            Printf("### track %d had to be killed due to crossing problems", fParticleV[i]);
            continue;
         }   
         nav->FindNextBoundaryAndStep(1.E30, kFALSE);
         fNextpathV[i]->InitFromNavigator(nav);
         fSnextV[i] += nav->GetStep();
      }
      if (fSnextV[i]>2.*gTolerance) fIzeroV[i] = 0;
   }
}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateTracks(GeantTrack_v &output)
{
// Propagate the ntracks with their selected steps. If a boundary is
// found in the way, the track is marked as hole. Nout must be initialized from outside.
//     trackin = array of <ntracks> input tracks
//     trackout = array of <nout> tracks propagated to physics processes
//     tracktodo = array of <ntodo> tracks propagated with a safe step or crossing 
//                 inside the same volume. These have to be propagated  again.
//     trackcross = array of <ncross> tracks that crossed the boundary. For these tracks
//                 the continuous processes have to be applied after propagation
   Int_t itr = 0;
   Int_t nsel = 0;
   Double_t step, snext, safety, c;
   const Double_t bmag = gPropagator->fBmag;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   Int_t tid = nav->GetThreadId();
   GeantThreadData *td = gPropagator->fThreadData[tid];
   // Remove dead tracks, propagate neutrals
   for (itr=0; itr<fNtracks; itr++) {
      if (fStatusV[itr] == kKilled) {
         MarkRemoved(itr);
         continue;
      }
      if (fChargeV[itr]==0 || bmag<1.E-10) {
      // Do straight propagation to physics process or boundary
         if (fFrombdrV[itr]) {
            *fPathV[itr] = *fNextPathV[itr];
            fStatus[itr] = kBoundary;
         } else {
            fStatus[itr] = kPhysics;
         }
         fPstepV[itr] -= fSnextV[itr];
         fSafetyV[itr] -= fSnextV[itr];
         if (fSafetyV[itr]<gTolerance) fSafetyV[itr] = 0;
         fXposV[itr] += fSnextV[itr]*fXdir[itr];
         fYposV[itr] += fSnextV[itr]*fYdir[itr];
         fZposV[itr] += fSnextV[itr]*fZdir[itr];
         fNsteps[itr]++;
         MarkRemoved(itr);
      }
   }
   // Compact remaining tracks and move the removed oned to the output container
   if (!fCompact) Compact(output);
   if (!fNtracks) return 0;

   // ONLY CHARGED TRACKS IN MAG. FIELD
   Double_t *steps = new Double_t[fNtracks];
   // Select tracks that undergo safe steps to physics process
   nsel = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (fPstep[itr]<fSafety[itr]) {
         Select(itr);
         fStatusV[itr] = kPhysics;
         nsel++;
      }   
   }
   // fPstep array is contiguous after the last Reshuffle operation
   if (nsel) {
      Reshuffle();
      memcpy(steps, fPstep, nsel*sizeof(Double_t));
      // This category of tracks can make the full step to the physics process
      PropagateInVolume(steps);
      // Move these tracks to the output container
      output.AddTracks(*this, 0, nsel-1);
      RemoveTracks(0, nsel-1);
   }
   // Check if we can propagate to boundary
   nsel = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      c = Curvature(itr);
      if (0.25*c*fSnextV[itr]<1E-6 && fSnextV[itr]<1E-3 && fSnextV[itr]<fStepV[itr]-1E-6) {
         Select[itr];
         nsel++;
      }
   } 
   if (nsel>0) {
      Reshuffle();
      for (Int_t itr=0; itr<nsel;      
         // Propagate with snext and check if we crossed
         //   backup track position and momentum
         if (track->izero>10) snext = 1.E-3;
         basket = track->PropagateInField(snext+10*gTolerance, kTRUE, trackin[itr]);
         if (snext<1.E-6) track->izero++;
         else track->izero = 0;
         gPropagator->fNsnextSteps++;
         if (!basket) {
            // track exiting
            if (nav->IsOutside()) {
               if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
               trackcross[ncross++] = trackin[itr];
               gPropagator->StopTrack(track);
               continue;
            }   
            // these tracks do not cross
//            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
            tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
            continue; // -> next track
         }
         // These tracks are reaching boundaries
         trackcross[ncross++] = trackin[itr];
//         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d pushed to boundary of %s", trackin[itr], basket->GetName());
//         basket = track->PropagateStraight(snext, trackin[itr]);
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
            gPropagator->fNsnextSteps++;
            if (!basket) {
               // track exiting geometry
               if (nav->IsOutside()) {
                  if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
                  trackcross[ncross++] = trackin[itr];
                  gPropagator->StopTrack(track);
                  continue;
               }   
               // these tracks do not cross
//               if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with snext=%19.15f", trackin[itr], snext);
               tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
               continue; // -> next track
            }
            // These tracks are reaching boundaries
            trackcross[ncross++] = trackin[itr];
            continue;
         }
         if (safety<1.E-3) track->izero++;
         // Propagate with safety without checking crossing
         basket = track->PropagateInField(safety, kFALSE, trackin[itr]);
      }   
      gPropagator->fNsafeSteps++;
      if (!basket) {
         // check if track exiting
         if (nav->IsOutside()) {
            if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d exiting geometry", trackin[itr]);
            trackcross[ncross++] = trackin[itr];
            gPropagator->StopTrack(track);
            continue;
         }   
         // these tracks do not cross
//         if (gPropagator->fUseDebug && (gPropagator->fDebugTrk== trackin[itr] || gPropagator->fDebugTrk<0)) Printf("   track %d propagated with safety=%19.15f", trackin[itr], safety);
         tracktodo[ntodo++] = trackin[itr]; // <- survives partial geometry step
         continue; // -> next track
      }
      // These tracks are reaching boundaries
      trackcross[ncross++] = trackin[itr];
   }   
   // Recompute snext and safety for todo tracks
   if (ntodo) ComputeTransportLength(ntodo, tracktodo);
}
