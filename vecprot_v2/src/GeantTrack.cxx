#include "globals.h"

#if USE_VECGEOM_NAVIGATOR == 1
#pragma message("Compiling against VecGeom")
#include "navigation/navigationstate.h"
#include "navigation/simple_navigator.h"
#include "volumes/placed_volume.h" // equivalent of TGeoNode
#include "base/vector3d.h"
#include "base/transformation_matrix.h"
#include "base/global.h"
#include "management/geo_manager.h"
typedef vecgeom::NavigationState VolumePath_t;
#include "TGeoBranchArray.h"
#include "TGeoNavigator.h"
#include "TGeoNode.h"

#else // TGeoNavigator as default
#pragma message("Compiling against TGeo")
#include "TGeoBranchArray.h"
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#include <iostream>
typedef TGeoBranchArray VolumePath_t;
#include "navigation/navigationstate.h"
#include "navigation/simple_navigator.h"
#include "volumes/placed_volume.h" // equivalent of TGeoNode
#include "base/vector3d.h"
#include "base/transformation_matrix.h"
#include "base/global.h"
#include "management/geo_manager.h"
#endif


#include "TGeoManager.h"
#include "GeantTrack.h"

#ifdef __STAT_DEBUG_TRK
#include "GeantTrackStat.h"
#endif

#include "GeantThreadData.h"
//#include "WorkloadManager.h"
#include "TGeoHelix.h"
#include "GeantScheduler.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h>
#else
#include "mm_malloc.h"
#endif
#include <cassert>

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
            fEindex(0),
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
            fEdep(0),
            fPstep(1.E20),
            fStep(0),
            fSnext(0),
            fSafety(0),
            fFrombdr(false),
            fPending(false),
            fPath(0),
            fNextpath(0)
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
            fEindex(0),
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
            fEdep(0),
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

#ifdef USE_VECGEOM_NAVIGATOR
      // At this time, the geometry is probably already
      // loaded and we could query it for the actual size
      fPath = new vecgeom::NavigationState( vecgeom::GeoManager::Instance().getMaxDepth() );
      fNextpath = new vecgeom::NavigationState( vecgeom::GeoManager::Instance().getMaxDepth() );
#else
      fPath = new TGeoBranchArray(30);
      fNextpath = new TGeoBranchArray(30);
#endif

}

//______________________________________________________________________________
GeantTrack::GeantTrack(const GeantTrack& other)
           :TObject(other),
            fEvent(other.fEvent),
            fEvslot(other.fEvslot),
            fParticle(other.fParticle),
            fPDG(other.fPDG),
            fG5code(other.fG5code),
            fEindex(other.fEindex),
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
            fEdep(other.fEdep),
            fPstep(other.fPstep),
            fStep(other.fStep),
            fSnext(other.fSnext),
            fSafety(other.fSafety),
            fFrombdr(other.fFrombdr),
            fPending(other.fPending),

            fPath(new VolumePath_t(*other.fPath)),
            fNextpath(new VolumePath_t(*other.fNextpath))
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
      fEindex = other.fEindex;
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
      fEdep = other.fEdep;
      fPstep = other.fPstep;
      fStep = other.fStep;
      fSnext = other.fSnext;
      fSafety = other.fSafety;
      fFrombdr = other.fFrombdr;
      fPending = other.fPending;
      fPath = new VolumePath_t(*other.fPath);
      fNextpath = new VolumePath_t(*other.fNextpath);
   }
   return *this;
}

//______________________________________________________________________________
GeantTrack::~GeantTrack()
{
// Destructor.
   delete fPath;
   delete fNextpath;
}

//______________________________________________________________________________
void GeantTrack::ReadFromVector(const GeantTrack_v &arr, Int_t i)
{
// Fill track from array
   fEvent = arr.fEventV[i];
   fEvslot = arr.fEvslotV[i];
   fParticle = arr.fParticleV[i];
   fPDG = arr.fPDGV[i];
   fG5code = arr.fG5codeV[i];
   fEindex = arr.fEindexV[i];
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
   fEdep = arr.fEdepV[i];
   fPstep = arr.fPstepV[i];
   fStep = arr.fStepV[i];
   fSnext = arr.fSnextV[i];
   fSafety = arr.fSafetyV[i];
   fFrombdr = arr.fFrombdrV[i];
   fPending = arr.fPendingV[i];

// if (fPath) *fPath = *arr.fPathV[i]; else fPath = new TGeoBranchArray(*arr.fPathV[i]);
// if (fNextpath) *fNextpath = *arr.fNextpathV[i]; else fNextpath = new TGeoBranchArray(*arr.fNextpathV[i]);
   if (fPath) *fPath = *arr.fPathV[i]; else fPath = new VolumePath_t(*arr.fPathV[i]);
   if (fNextpath) *fNextpath = *arr.fNextpathV[i]; else fNextpath = new VolumePath_t(*arr.fNextpathV[i]);
}

//______________________________________________________________________________
void GeantTrack::Clear(Option_t *)
{
// Resets track content.
   fEvent = -1;
   fEvslot = -1;
   fParticle = -1;
   fPDG = 0;
   fG5code = 0;
   fEindex = 0;
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
   fEdep = 0;
   fPstep = 1.E20;
   fStep = 0.;
   fSnext = 0.;
   fSafety = 0.;
   fFrombdr = false;
   fPending = false;
}

//______________________________________________________________________________
Double_t GeantTrack::Curvature() const
{
// Curvature
   if (fCharge==0) return 0.;
   return TMath::Abs(kB2C*gPropagator->fBmag/Pt());
}


//______________________________________________________________________________
void GeantTrack::SetPath(VolumePath_t const * const path)
{
   *fPath = *path;
}

//______________________________________________________________________________
void GeantTrack::SetNextPath(VolumePath_t const * const path)
{
// Set next path.
   *fNextpath = *path;
}


//______________________________________________________________________________
void GeantTrack::Print(Int_t) const {
   TString spath;
//   if (path) path->GetPath(spath);
   Printf("=== Track %d (ev=%d): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Dir:(%f,%f,%f) P:%g E:%g snext=%g safety=%g nsteps=%d",
           fParticle,fEvent, fProcess,fPstep,fCharge,fXpos,fYpos,fZpos,fXdir,fYdir,fZdir,P(),fE,fSnext,fSafety,fNsteps);
   Printf("G5Code %d, Mass %lf\n", fG5code, fMass);
}


ClassImp(GeantTrack_v)

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v()
             :fNtracks(0),fMaxtracks(0),fNselected(0),fHoles(),fSelected(),fCompact(true),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fEindexV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fEdepV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
   // Printf("GeantTrack_v default %p",this);
    // Dummy ctor.
#ifdef __STAT_DEBUG_TRK
   fStat.InitArrays(gPropagator->fNevents);
#endif
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(Int_t size)
             :fNtracks(0),fMaxtracks(0),fNselected(0),fHoles(size),fSelected(size),fCompact(true),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fEindexV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fEdepV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
    //Printf("GeantTrack_v default2 %p",this);
    // Constructor with maximum capacity.
#ifdef __STAT_DEBUG_TRK
   fStat.InitArrays(gPropagator->fNevents);
#endif
   Resize(size);
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(const GeantTrack_v &track_v)
             :fNtracks(track_v.fNtracks),fMaxtracks(track_v.fMaxtracks),fNselected(track_v.fNselected),fHoles(track_v.fHoles),
              fSelected(track_v.fSelected),fCompact(track_v.fCompact),fBuf(0),fEventV(0),fEvslotV(0),fParticleV(0),
              fPDGV(0),fG5codeV(0),fEindexV(0),fChargeV(0),fProcessV(0),fIzeroV(0),fNstepsV(0),
              fSpeciesV(0),fStatusV(0),fMassV(0),fXposV(0),fYposV(0),fZposV(0),
              fXdirV(0),fYdirV(0),fZdirV(0),fPV(0),fEV(0),fEdepV(0),fPstepV(0),fStepV(0),
              fSnextV(0),fSafetyV(0),fFrombdrV(0),fPendingV(0),fPathV(0),fNextpathV(0)
{
// Copy constructor
   // Printf("GeantTrack_v copy const %p",this);
    #ifdef __STAT_DEBUG_TRK
   fStat.InitArrays(gPropagator->fNevents);
#endif
   Int_t size = track_v.fMaxtracks;
   fBuf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
   memset(fBuf, 0, size*sizeof(GeantTrack));
   memcpy(fBuf, track_v.fBuf, fNtracks*sizeof(GeantTrack));
   AssignInBuffer(fBuf, size);
}

//______________________________________________________________________________
GeantTrack_v &GeantTrack_v::operator=(const GeantTrack_v &track_v)
{
    //Printf("GeantTrack_v assignment %p",this);
    // Assignment operator
   if (&track_v != this) {
      Int_t size = track_v.fMaxtracks;
      if (fMaxtracks<size) {
         _mm_free(fBuf);
         fBuf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
      }
      fNtracks = track_v.fNtracks;
      fMaxtracks = size;
      fNselected = track_v.fNselected;
      fHoles = track_v.fHoles;
      fSelected = track_v.fSelected;
      fCompact = track_v.fCompact;
      memcpy(fBuf, track_v.fBuf, size*sizeof(GeantTrack));
      AssignInBuffer(fBuf, size);
#ifdef __STAT_DEBUG_TRK
   fStat.InitArrays(gPropagator->fNevents);
#endif
   }
   return *this;
}

//______________________________________________________________________________
GeantTrack_v::~GeantTrack_v()
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
   fEindexV = (Int_t*)buf;
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
   fEdepV = (Double_t*)buf;
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

//   fPathV = (TGeoBranchArray**)buf;
   //buf += size*sizeof(TGeoBranchArray*);
   //fNextpathV = (TGeoBranchArray**)buf;
   //buf += size*sizeof(TGeoBranchArray*);
   fPathV = (VolumePath_t**)buf;
   buf += size*sizeof(VolumePath_t*);
   fNextpathV = (VolumePath_t**)buf;
   buf += size*sizeof(VolumePath_t*);
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
   memcpy(buf, fEventV, size_int);
   fEventV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fEvslotV, size_int);
   fEvslotV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fParticleV, size_int);
   fParticleV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fPDGV, size_int);
   fPDGV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fG5codeV, size_int);
   fG5codeV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fEindexV, size_int);
   fEindexV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fChargeV, size_int);
   fChargeV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fProcessV, size_int);
   fProcessV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fIzeroV, size_int);
   fIzeroV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fNstepsV, size_int);
   fNstepsV = (Int_t*)buf;
   buf += size_intn;
   memcpy(buf, fSpeciesV, fNtracks*sizeof(Species_t));
   fSpeciesV = (Species_t*)buf;
   buf += size*sizeof(Species_t);
   memcpy(buf, fStatusV, fNtracks*sizeof(TrackStatus_t));
   fStatusV = (TrackStatus_t*)buf;
   buf += size*sizeof(TrackStatus_t);
   memcpy(buf, fMassV, size_double);
   fMassV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fXposV, size_double);
   fXposV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fYposV, size_double);
   fYposV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fZposV, size_double);
   fZposV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fXdirV, size_double);
   fXdirV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fYdirV, size_double);
   fYdirV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fZdirV, size_double);
   fZdirV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fPV, size_double);
   fPV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fEV, size_double);
   fEV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fEdepV, size_double);
   fEdepV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fPstepV, size_double);
   fPstepV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fStepV, size_double);
   fStepV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fSnextV, size_double);
   fSnextV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fSafetyV, size_double);
   fSafetyV = (Double_t*)buf;
   buf += size_doublen;
   memcpy(buf, fFrombdrV, fNtracks*sizeof(Bool_t));
   fFrombdrV = (Bool_t*)buf;
   buf += size*sizeof(Bool_t);
   memcpy(buf, fPendingV, fNtracks*sizeof(Bool_t));
   fPendingV = (Bool_t*)buf;
   buf += size*sizeof(Bool_t);
   // Eventhough the fPath are pointers, this is fine as CopyToBuffer
   // is used to resize the underlying arrays and the previous ones
   // are just dropped.
   // However, because we need fPathV to be initalized we copy
   // past the end of the active part of the array.

//   memcpy(buf, fPathV, fMaxtracks*sizeof(TGeoBranchArray*));
   //fPathV = (TGeoBranchArray**)buf;
  // buf += size*sizeof(TGeoBranchArray*);
  // memcpy(buf, fNextpathV, fMaxtracks*sizeof(TGeoBranchArray*));
  // fNextpathV = (TGeoBranchArray**)buf;
  // buf += size*sizeof(TGeoBranchArray*);

   memcpy(buf, fPathV, fMaxtracks*sizeof(VolumePath_t*));
   fPathV = (VolumePath_t**)buf;
   buf += size*sizeof(VolumePath_t*);
   memcpy(buf, fNextpathV, fMaxtracks*sizeof(VolumePath_t*));
   fNextpathV = (VolumePath_t**)buf;
   buf += size*sizeof(VolumePath_t*);
}

//______________________________________________________________________________
Bool_t GeantTrack_v::IsSame(const GeantTrack_v &tr1, Int_t i1, const GeantTrack_v &tr2, Int_t i2)
{
// Compare two tracks.
   Long64_t chk1, chk2;
   chk1 = tr1.fEventV[i1] + tr1.fEvslotV[i1] + tr1.fParticleV[i1] + tr1.fPDGV[i1]
        + tr1.fG5codeV[i1] + tr1.fEindexV[i1] + tr1.fChargeV[i1] + tr1.fProcessV[i1] + tr1.fIzeroV[i1]
        + tr1.fNstepsV[i1] + (Long64_t)tr1.fSpeciesV[i1] + (Long64_t)tr1.fStatusV[i1];
   chk2 = tr2.fEventV[i2] + tr2.fEvslotV[i2] + tr2.fParticleV[i2] + tr2.fPDGV[i2]
        + tr2.fG5codeV[i2] + tr2.fEindexV[i2] + tr2.fChargeV[i2] + tr2.fProcessV[i2] + tr2.fIzeroV[i2]
        + tr2.fNstepsV[i2] + (Long64_t)tr2.fSpeciesV[i2] + (Long64_t)tr2.fStatusV[i2];
   if (chk1 != chk2){
    Printf("difference in integer properties\n");
    return false;
   }
   Double_t dchk1, dchk2;
   dchk1 = (Long64_t)tr1.fMassV[i1] + tr1.fXposV[i1] + tr1.fYposV[i1] + tr1.fZposV[i1]
        + tr1.fXdirV[i1] + tr1.fYdirV[i1] + tr1.fZdirV[i1] + tr1.fPV[i1] + tr1.fEdepV[i1]
        + tr1.fEV[i1] + tr1.fPstepV[i1] + tr1.fStepV[i1] + tr1.fSnextV[i1] + tr1.fSafetyV[i1];
   dchk2 = (Long64_t)tr2.fMassV[i2] + tr2.fXposV[i2] + tr2.fYposV[i2] + tr2.fZposV[i2]
        + tr2.fXdirV[i2] + tr2.fYdirV[i2] + tr2.fZdirV[i2] + tr2.fPV[i2] + tr2.fEdepV[i2]
        + tr2.fEV[i2] + tr2.fPstepV[i2] + tr2.fStepV[i2] + tr2.fSnextV[i2] + tr2.fSafetyV[i2];
   if (!TMath::AreEqualAbs(dchk1, dchk2, 1.E-10))
   {
       Printf("difference in double properties; delta %lf \n", dchk1-dchk2);
       Printf("Mass %lf %lf", tr1.fMassV[i1], tr2.fMassV[i1] );
       Printf("fXPosV %lf %lf", tr1.fXposV[i1], tr2.fXposV[i1]  );
       Printf("fYPosV %lf %lf", tr1.fYposV[i1], tr2.fYposV[i1] );
       Printf("fZPosV %lf %lf", tr1.fZposV[i1], tr2.fZposV[i1] );
       Printf("fXDirV %lf %lf", tr1.fXdirV[i1], tr2.fXdirV[i1] );
       Printf("fYDirV %lf %lf", tr1.fYdirV[i1], tr2.fYdirV[i1] );
       Printf("fZdirV %lf %lf", tr1.fZdirV[i1], tr2.fZdirV[i1] );
       Printf("fPV %lf %lf", tr1.fPV[i1], tr2.fPV[i1] );
       Printf("fEdepV %lf %lf", tr1.fEdepV[i1], tr2.fEdepV[i1] );
       Printf("fEV %lf %lf", tr1.fEV[i1], tr2.fEV[i1] );
       Printf("fPStepV %lf %lf", tr1.fPstepV[i1], tr2.fPstepV[i1] );
       Printf("StepV %lf %lf", tr1.fStepV[i1], tr2.fStepV[i1] );
       Printf("Snext %lf %lf", tr1.fSnextV[i1],tr2.fSnextV[i1] );
       Printf("Safety %lf %lf", tr1.fSafetyV[i1], tr2.fSafetyV[i1] );
       return false;
   }
   if (tr1.fPendingV[i1] != tr2.fPendingV[i2])
   {
       Printf("difference in scheduling properties\n");
       return false;
   }
   return true;
}

//______________________________________________________________________________
void GeantTrack_v::CheckTracks()
{
 //  for (Int_t i=0; i<fMaxtracks; ++i)
 //     if (fNextpathV[i]) fNextpathV[i]->SetClient(gPropagator);
}

//______________________________________________________________________________
void GeantTrack_v::Resize(Int_t newsize)
{
   // Printf("Resizing %p", this);
    // Resize the container.
   Int_t size = round_up_align(newsize);
   if (size<fNtracks) {
      printf("Error: Cannot resize to less than current track content\n");
      return;
   }
   fHoles.ResetBitNumber(size-1);
   fSelected.ResetBitNumber(size-1);
   char* buf = (char*)_mm_malloc(size*sizeof(GeantTrack), ALIGN_PADDING);
   memset(buf, 0, size*sizeof(GeantTrack));
   if (!fBuf) {
      // All arrays are contiguous in a single buffer and aligned with the
      // same padding ALIGN_PADDING
      fBuf = buf;
      fMaxtracks = size;
      AssignInBuffer(buf, size);

//      memset( fPathV, 0, size * sizeof(TGeoBranchArray*) );
      //memset( fNextpathV, 0, size * sizeof(TGeoBranchArray*) );
      memset( fPathV, 0, size * sizeof(VolumePath_t*) );
      memset( fNextpathV, 0, size * sizeof(VolumePath_t*) );
   } else {
      // Resize container
      CopyToBuffer(buf, size);
      _mm_free(fBuf);
      fBuf = buf;
      UInt_t increase = size - fMaxtracks;
//      memset( fPathV+fMaxtracks, 0, increase * sizeof(TGeoBranchArray*) );
//      memset( fNextpathV+fMaxtracks, 0, increase * sizeof(TGeoBranchArray*) );
      memset( fPathV+fMaxtracks, 0, increase * sizeof(VolumePath_t*) );
      memset( fNextpathV+fMaxtracks, 0, increase * sizeof(VolumePath_t*) );
      fMaxtracks = size;
   }
}

//______________________________________________________________________________
Int_t GeantTrack_v::AddTrack(const GeantTrack &track)
{
#ifdef VERBOSE
   Printf("AddTrack on %p from %p",this, &track);
   track.Print();
#endif
    // Add new track to the array. If addition is done on top of non-compact array,
   // the track will be inserted without updating the number of tracks.
   // Returns the location where the track was added.
   Int_t itrack = fNtracks;
   if (!fCompact) itrack = fHoles.FirstSetBit();
   if (itrack==fMaxtracks) Resize(2*fMaxtracks);
   fHoles.ResetBitNumber(itrack);
   fSelected.ResetBitNumber(itrack);
   fEventV    [itrack] = track.fEvent;
   fEvslotV   [itrack] = track.fEvslot;
   fParticleV [itrack] = track.fParticle;
   fPDGV      [itrack] = track.fPDG;
   fG5codeV   [itrack] = track.fG5code;
   fEindexV   [itrack] = track.fEindex;
   fChargeV   [itrack] = track.fCharge;
   fProcessV  [itrack] = track.fProcess;
   fIzeroV    [itrack] = track.fIzero;
   fNstepsV   [itrack] = track.fNsteps;
   fSpeciesV  [itrack] = track.fSpecies;
   fStatusV   [itrack] = track.fStatus;
   fMassV     [itrack] = track.fMass;
   fXposV     [itrack] = track.fXpos;
   fYposV     [itrack] = track.fYpos;
   fZposV     [itrack] = track.fZpos;
   fXdirV     [itrack] = track.fXdir;
   fYdirV     [itrack] = track.fYdir;
   fZdirV     [itrack] = track.fZdir;
   fPV        [itrack] = track.fP;
   fEV        [itrack] = track.fE;
   fEdepV     [itrack] = track.fEdep;
   fPstepV    [itrack] = track.fPstep;
   fStepV     [itrack] = track.fStep;
   fSnextV    [itrack] = track.fSnext;
   fSafetyV   [itrack] = track.fSafety;
   fFrombdrV  [itrack] = track.fFrombdr;
   fPendingV  [itrack] = track.fPending;

//   if (fPathV[itrack]) *fPathV[itrack] = *track.fPath; else fPathV[itrack] = new TGeoBranchArray(*track.fPath);
//   if (fNextpathV[itrack]) *fNextpathV[itrack] = *track.fNextpath; else fNextpathV[itrack] = new TGeoBranchArray(*track.fNextpath);

   if (fPathV[itrack]) *fPathV[itrack] = *track.fPath; else fPathV[itrack] = new VolumePath_t(*track.fPath);
   if (fNextpathV[itrack]) *fNextpathV[itrack] = *track.fNextpath; else fNextpathV[itrack] = new VolumePath_t(*track.fNextpath);

   if (itrack==fNtracks) fNtracks++;
#ifdef __STAT_DEBUG_TRK
   fStat.fNtracks[fEvslotV[itrack]]++;
#endif

   return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::GetTrack(Int_t i, GeantTrack &track) const {
    Printf("GetTrack %p",this);
    // Extract a single track from array.
   track.ReadFromVector(*this, i);
}


//______________________________________________________________________________
Int_t GeantTrack_v::AddTrack(const GeantTrack_v &arr, Int_t i)
{
   //Printf("AddTrack %d from %p to %p",i, &arr, this);
   // Add track from different array
   // If addition is done on top of non-compact array,
   // the track will be inserted without updating the number of tracks.
   // Returns the location where the track was added.
#ifdef VERBOSE
	arr.PrintTrack( i );
#endif
	if( ( arr.fG5codeV[i]==0 && arr.fMassV[i]==0 ) )
     {
          arr.PrintTrack( i );
          assert( ! ( arr.fG5codeV[i]==0 && arr.fMassV[i]==0 ) );
      }
   Int_t itrack = fNtracks;
   if (!fCompact) itrack = fHoles.FirstSetBit();
   if (itrack==fMaxtracks) Resize(2*fMaxtracks);
   fHoles.ResetBitNumber(itrack);

   fEventV    [fNtracks] = arr.fEventV    [i];
   fEvslotV   [fNtracks] = arr.fEvslotV   [i];
   fParticleV [fNtracks] = arr.fParticleV [i];
   fPDGV      [fNtracks] = arr.fPDGV      [i];
   fG5codeV   [fNtracks] = arr.fG5codeV   [i];
   fEindexV   [fNtracks] = arr.fEindexV   [i];
   fChargeV   [fNtracks] = arr.fChargeV   [i];
   fProcessV  [fNtracks] = arr.fProcessV  [i];
   fIzeroV    [fNtracks] = arr.fIzeroV    [i];
   fNstepsV   [fNtracks] = arr.fNstepsV   [i];
   fSpeciesV  [fNtracks] = arr.fSpeciesV  [i];
   fStatusV   [fNtracks] = arr.fStatusV   [i];
   fMassV     [fNtracks] = arr.fMassV     [i];
   fXposV     [fNtracks] = arr.fXposV     [i];
   fYposV     [fNtracks] = arr.fYposV     [i];
   fZposV     [fNtracks] = arr.fZposV     [i];
   fXdirV     [fNtracks] = arr.fXdirV     [i];
   fYdirV     [fNtracks] = arr.fYdirV     [i];
   fZdirV     [fNtracks] = arr.fZdirV     [i];
   fPV        [fNtracks] = arr.fPV        [i];
   fEV        [fNtracks] = arr.fEV        [i];
   fEdepV     [fNtracks] = arr.fEdepV     [i];
   fPstepV    [fNtracks] = arr.fPstepV    [i];
   fStepV     [fNtracks] = arr.fStepV     [i];
   fSnextV    [fNtracks] = arr.fSnextV    [i];
   fSafetyV   [fNtracks] = arr.fSafetyV   [i];
   fFrombdrV  [fNtracks] = arr.fFrombdrV  [i];
   fPendingV  [fNtracks] = arr.fPendingV  [i];
   if (fPathV[fNtracks]) *fPathV[fNtracks] = *arr.fPathV[i];
   else
//	   fPathV[fNtracks] = new TGeoBranchArray(*arr.fPathV[i]);
	   fPathV[fNtracks] = new VolumePath_t(*arr.fPathV[i]);
   if (fNextpathV[fNtracks]) *fNextpathV[fNtracks] = *arr.fNextpathV[i];
   else
	   // fNextpathV[fNtracks] = new TGeoBranchArray(*arr.fNextpathV[i]);
	   fNextpathV[fNtracks] = new VolumePath_t(*arr.fNextpathV[i]);
   fSelected.ResetBitNumber(fNtracks);
   fHoles.ResetBitNumber(fNtracks);

   if (itrack==fNtracks) fNtracks++;
   if (TMath::IsNaN(fXdirV[fNtracks-1]) || !IsSame(arr,i, *this, fNtracks-1)) {
      Printf("Error: AddTrack(GeantTrack const &): Different tracks");
      Printf("This happened for index %d\n",i);
      Printf("IsNan check %d\n", TMath::IsNaN(fXdirV[fNtracks-1]));
      Printf("difference in checksum %d\n", !IsSame(arr,i, *this, fNtracks-1));
   }
#ifdef __STAT_DEBUG_TRK
   fStat.fNtracks[arr.fEvslotV[i]]++;
#endif

   return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::AddTracks(const GeantTrack_v &arr, Int_t istart, Int_t iend)
{
   // Printf("AddTrack %d to %d from %p to %p",istart, iend, &arr, this);
    if( ( arr.fG5codeV[istart]==0 && arr.fMassV[istart]==0 ) )
    {
        arr.PrintTrack( istart );
        assert( ! ( arr.fG5codeV[istart]==0 && arr.fMassV[istart]==0 ) );
    }
    // Add tracks from different array
#ifdef __STAT_DEBUG_TRK
   for (Int_t i=istart; i<=iend; i++) fStat.fNtracks[arr.fEvslotV[i]]++;
#endif
   Int_t ncpy = iend-istart+1;
   if (fNtracks+ncpy>=fMaxtracks) Resize(TMath::Max(2*fMaxtracks, fNtracks+ncpy));
   memcpy(&fEventV    [fNtracks], &arr.fEventV    [istart], ncpy*sizeof(Int_t));
   memcpy(&fEvslotV   [fNtracks], &arr.fEvslotV   [istart], ncpy*sizeof(Int_t));
   memcpy(&fParticleV [fNtracks], &arr.fParticleV [istart], ncpy*sizeof(Int_t));
   memcpy(&fPDGV      [fNtracks], &arr.fPDGV      [istart], ncpy*sizeof(Int_t));
   memcpy(&fG5codeV   [fNtracks], &arr.fG5codeV   [istart], ncpy*sizeof(Int_t));
   memcpy(&fEindexV   [fNtracks], &arr.fEindexV   [istart], ncpy*sizeof(Int_t));
   memcpy(&fChargeV   [fNtracks], &arr.fChargeV   [istart], ncpy*sizeof(Int_t));
   memcpy(&fProcessV  [fNtracks], &arr.fProcessV  [istart], ncpy*sizeof(Int_t));
   memcpy(&fIzeroV    [fNtracks], &arr.fIzeroV    [istart], ncpy*sizeof(Int_t));
   memcpy(&fNstepsV   [fNtracks], &arr.fNstepsV   [istart], ncpy*sizeof(Int_t));
   memcpy(&fSpeciesV  [fNtracks], &arr.fSpeciesV  [istart], ncpy*sizeof(Species_t));
   memcpy(&fStatusV   [fNtracks], &arr.fStatusV   [istart], ncpy*sizeof(TrackStatus_t));
   memcpy(&fMassV     [fNtracks], &arr.fMassV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fXposV     [fNtracks], &arr.fXposV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fYposV     [fNtracks], &arr.fYposV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fZposV     [fNtracks], &arr.fZposV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fXdirV     [fNtracks], &arr.fXdirV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fYdirV     [fNtracks], &arr.fYdirV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fZdirV     [fNtracks], &arr.fZdirV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fPV        [fNtracks], &arr.fPV        [istart], ncpy*sizeof(Double_t));
   memcpy(&fEV        [fNtracks], &arr.fEV        [istart], ncpy*sizeof(Double_t));
   memcpy(&fEdepV     [fNtracks], &arr.fEdepV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fPstepV    [fNtracks], &arr.fPstepV    [istart], ncpy*sizeof(Double_t));
   memcpy(&fStepV     [fNtracks], &arr.fStepV     [istart], ncpy*sizeof(Double_t));
   memcpy(&fSnextV    [fNtracks], &arr.fSnextV    [istart], ncpy*sizeof(Double_t));
   memcpy(&fSafetyV   [fNtracks], &arr.fSafetyV   [istart], ncpy*sizeof(Double_t));
   memcpy(&fFrombdrV  [fNtracks], &arr.fFrombdrV  [istart], ncpy*sizeof(Bool_t));
   memcpy(&fPendingV  [fNtracks], &arr.fPendingV  [istart], ncpy*sizeof(Bool_t));
   for(Int_t i = fNtracks, j = istart; i < (fNtracks+ncpy) ; ++i) {
      // The following was wrong ... because the fPath are pointers.
      // memcpy(&fPathV[fNtracks], &arr.fPathV[istart], ncpy*sizeof(TGeoBranchArray*));
      // memcpy(&fNextpathV[fNtracks], &arr.fNextpathV[istart], ncpy*sizeof(TGeoBranchArray*));
      if (fPathV[i]) *fPathV[i] = *arr.fPathV[j];
      //else fPathV[i]  = new TGeoBranchArray(*arr.fPathV[j]);
      else fPathV[i]  = new VolumePath_t(*arr.fPathV[j]);
      if (fNextpathV[i]) *fNextpathV[i] = *arr.fNextpathV[j];
      //else fNextpathV[i] = new TGeoBranchArray(*arr.fNextpathV[j]);
      else fNextpathV[i] = new VolumePath_t(*arr.fNextpathV[j]);
   }
   fSelected.ResetBitNumber(fNtracks+ncpy-1);
   fHoles.ResetBitNumber(fNtracks+ncpy-1);
   for (Int_t i=istart; i<=iend; i++){
          if (!IsSame(arr,i, *this, fNtracks+i-istart)) Printf("Error: AddTracks(arr): Different tracks");
   }
   fNtracks += ncpy;
}

//______________________________________________________________________________
void GeantTrack_v::SwapTracks(Int_t i, Int_t j)
{
//    Printf("swapping tracks\n");
    // Swap two tracks in the container
   Double_t tdbl;
   Int_t    tint;
   Bool_t   tbool;
   VolumePath_t *tptr;
   tint = fEventV    [i]; fEventV    [i] = fEventV    [j]; fEventV    [j] = tint;
   tint = fEvslotV   [i]; fEvslotV   [i] = fEvslotV   [j]; fEvslotV   [j] = tint;
   tint = fParticleV [i]; fParticleV [i] = fParticleV [j]; fParticleV [j] = tint;
   tint = fPDGV      [i]; fPDGV      [i] = fPDGV      [j]; fPDGV      [j] = tint;
   tint = fG5codeV   [i]; fG5codeV   [i] = fG5codeV   [j]; fG5codeV   [j] = tint;
   tint = fEindexV   [i]; fEindexV   [i] = fEindexV   [j]; fEindexV   [j] = tint;
   tint = fChargeV   [i]; fChargeV   [i] = fChargeV   [j]; fChargeV   [j] = tint;
   tint = fProcessV  [i]; fProcessV  [i] = fProcessV  [j]; fProcessV  [j] = tint;
   tint = fIzeroV    [i]; fIzeroV    [i] = fIzeroV    [j]; fIzeroV    [j] = tint;
   tint = fNstepsV   [i]; fNstepsV   [i] = fNstepsV   [j]; fNstepsV   [j] = tint;
   Species_t tspec =
          fSpeciesV  [i]; fSpeciesV  [i] = fSpeciesV  [j]; fSpeciesV  [j] = tspec;
   TrackStatus_t  tstat =
          fStatusV   [i]; fStatusV   [i] = fStatusV   [j]; fStatusV   [j] = tstat;
   tdbl = fMassV     [i]; fMassV     [i] = fMassV     [j]; fMassV     [j] = tdbl;
   tdbl = fXposV     [i]; fXposV     [i] = fXposV     [j]; fXposV     [j] = tdbl;
   tdbl = fYposV     [i]; fYposV     [i] = fYposV     [j]; fYposV     [j] = tdbl;
   tdbl = fZposV     [i]; fZposV     [i] = fZposV     [j]; fZposV     [j] = tdbl;
   tdbl = fXdirV     [i]; fXdirV     [i] = fXdirV     [j]; fXdirV     [j] = tdbl;
   tdbl = fYdirV     [i]; fYdirV     [i] = fYdirV     [j]; fYdirV     [j] = tdbl;
   tdbl = fZdirV     [i]; fZdirV     [i] = fZdirV     [j]; fZdirV     [j] = tdbl;
   tdbl = fPV        [i]; fPV        [i] = fPV        [j]; fPV        [j] = tdbl;
   tdbl = fEV        [i]; fEV        [i] = fEV        [j]; fEV        [j] = tdbl;
   tdbl = fEdepV     [i]; fEdepV     [i] = fEdepV     [j]; fEdepV     [j] = tdbl;
   tdbl = fPstepV    [i]; fPstepV    [i] = fPstepV    [j]; fPstepV    [j] = tdbl;
   tdbl = fStepV     [i]; fStepV     [i] = fStepV     [j]; fStepV     [j] = tdbl;
   tdbl = fSnextV    [i]; fSnextV    [i] = fSnextV    [j]; fSnextV    [j] = tdbl;
   tdbl = fSafetyV   [i]; fSafetyV   [i] = fSafetyV   [j]; fSafetyV   [j] = tdbl;
   tbool = fFrombdrV [i]; fFrombdrV  [i] = fFrombdrV  [j]; fFrombdrV  [j] = tbool;
   tbool = fPendingV [i]; fPendingV  [i] = fPendingV  [j]; fPendingV  [j] = tbool;
   tptr = fPathV     [i]; fPathV     [i] = fPathV     [j]; fPathV     [j] = tptr;
   tptr = fNextpathV [i]; fNextpathV [i] = fNextpathV [j]; fNextpathV [j] = tptr;
   Bool_t sel = fSelected.TestBitNumber(j);
   fSelected.SetBitNumber(j, fSelected.TestBitNumber(i));
   fSelected.SetBitNumber(i, sel);
}

//______________________________________________________________________________
void GeantTrack_v::ReplaceTrack(Int_t i, Int_t j)
{
  //  Printf("Replacing tracks %d by %d on %p",i,j,this);
    // Replace content of track i with the one of track j
   fEventV    [i] = fEventV    [j];
   fEvslotV   [i] = fEvslotV   [j];
   fParticleV [i] = fParticleV [j];
   fPDGV      [i] = fPDGV      [j];
   fG5codeV   [i] = fG5codeV   [j];
   fEindexV   [i] = fEindexV   [j];
   fChargeV   [i] = fChargeV   [j];
   fProcessV  [i] = fProcessV  [j];
   fIzeroV    [i] = fIzeroV    [j];
   fNstepsV   [i] = fNstepsV   [j];
   fSpeciesV  [i] = fSpeciesV  [j];
   fStatusV   [i] = fStatusV   [j];
   fMassV     [i] = fMassV     [j];
   fXposV     [i] = fXposV     [j];
   fYposV     [i] = fYposV     [j];
   fZposV     [i] = fZposV     [j];
   fXdirV     [i] = fXdirV     [j];
   fYdirV     [i] = fYdirV     [j];
   fZdirV     [i] = fZdirV     [j];
   fPV        [i] = fPV        [j];
   fEV        [i] = fEV        [j];
   fEdepV     [i] = fEdepV     [j];
   fPstepV    [i] = fPstepV    [j];
   fStepV     [i] = fStepV     [j];
   fSnextV    [i] = fSnextV    [j];
   fSafetyV   [i] = fSafetyV   [j];
   fFrombdrV  [i] = fFrombdrV  [j];
   fPendingV  [i] = fPendingV  [j];
   if (fPathV[i]) *fPathV[i] = *fPathV[j];
   //else fPathV[i] = new TGeoBranchArray(*fPathV[j]);
   else fPathV[i] = new VolumePath_t(*fPathV[j]);
   if (fNextpathV[i]) *fNextpathV[i] = *fNextpathV[j];
   //else fNextpathV[i] = new TGeoBranchArray(*fNextpathV[j]);
   else fNextpathV[i] = new VolumePath_t(*fNextpathV[j]);
   fSelected.SetBitNumber(i, fSelected.TestBitNumber(j));
}

//______________________________________________________________________________
void GeantTrack_v::DeleteTrack(Int_t itr)
{
    //Printf("Deleting tracks on %p",this);
    // Delete branch arrays for this track. The track should not have a copy, this has
// to be called after a killed track is removed by the scheduler.
   delete fPathV[itr];     fPathV[itr] = 0;
   delete fNextpathV[itr]; fNextpathV[itr] = 0;
   // MarkRemoved(itr);
}

//______________________________________________________________________________
void GeantTrack_v::RemoveTracks(Int_t from, Int_t to)
{
#ifdef VERBOSE
	Printf("Removing tracks on %p from %d to %d",this, from, to);
#endif
	// Remove tracks from the container
#ifdef __STAT_DEBUG_TRK
   for (Int_t i=from; i<=to; i++) fStat.fNtracks[fEvslotV[i]]--;
#endif
   Int_t ncpy = fNtracks-to;
   memmove(&fEventV    [from], &fEventV    [to+1], ncpy*sizeof(Int_t));
   memmove(&fEvslotV   [from], &fEvslotV   [to+1], ncpy*sizeof(Int_t));
   memmove(&fParticleV [from], &fParticleV [to+1], ncpy*sizeof(Int_t));
   memmove(&fPDGV      [from], &fPDGV      [to+1], ncpy*sizeof(Int_t));
   memmove(&fG5codeV   [from], &fG5codeV   [to+1], ncpy*sizeof(Int_t));
   memmove(&fEindexV   [from], &fEindexV   [to+1], ncpy*sizeof(Int_t));
   memmove(&fChargeV   [from], &fChargeV   [to+1], ncpy*sizeof(Int_t));
   memmove(&fProcessV  [from], &fProcessV  [to+1], ncpy*sizeof(Int_t));
   memmove(&fIzeroV    [from], &fIzeroV    [to+1], ncpy*sizeof(Int_t));
   memmove(&fNstepsV   [from], &fNstepsV   [to+1], ncpy*sizeof(Int_t));
   memmove(&fSpeciesV  [from], &fSpeciesV  [to+1], ncpy*sizeof(Species_t));
   memmove(&fStatusV   [from], &fStatusV   [to+1], ncpy*sizeof(TrackStatus_t));
   memmove(&fMassV     [from], &fMassV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fXposV     [from], &fXposV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fYposV     [from], &fYposV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fZposV     [from], &fZposV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fXdirV     [from], &fXdirV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fYdirV     [from], &fYdirV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fZdirV     [from], &fZdirV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fPV        [from], &fPV        [to+1], ncpy*sizeof(Double_t));
   memmove(&fEV        [from], &fEV        [to+1], ncpy*sizeof(Double_t));
   memmove(&fEdepV     [from], &fEdepV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fPstepV    [from], &fPstepV    [to+1], ncpy*sizeof(Double_t));
   memmove(&fStepV     [from], &fStepV     [to+1], ncpy*sizeof(Double_t));
   memmove(&fSnextV    [from], &fSnextV    [to+1], ncpy*sizeof(Double_t));
   memmove(&fSafetyV   [from], &fSafetyV   [to+1], ncpy*sizeof(Double_t));
   memmove(&fFrombdrV  [from], &fFrombdrV  [to+1], ncpy*sizeof(Bool_t));
   memmove(&fPendingV  [from], &fPendingV  [to+1], ncpy*sizeof(Bool_t));
   for (Int_t i = from, j = to+1, k = 0; k < ncpy; ++i,++j,++k) {
      // This was wrong, we must delete the overwritten one and zero the 'moved' part,
      // or we need to swap them.  (This code has memory leaks and double use/delete)
      // memmove(&fPathV[from], &fPathV[to+1], ncpy*sizeof(TGeoBranchArray*));
      // memmove(&fNextpathV[from], &fNextpathV[to+1], ncpy*sizeof(TGeoBranchArray*));
      //TGeoBranchArray *tptr;
      VolumePath_t *tptr;
	  tptr = fPathV[i]; fPathV[i] = fPathV[j]; fPathV[j] = tptr;
      tptr = fNextpathV[i]; fNextpathV[i] = fNextpathV[j]; fNextpathV[j] = tptr;
   }
   fNtracks -= to-from+1;
   fSelected.ResetAllBits();
   fNselected = 0;
}

//______________________________________________________________________________
Int_t GeantTrack_v::Compact(GeantTrack_v *moveto)
{
    //Printf("Compacting tracks on %p",this);
    // Compact the holes in the array. Return number of active elements. This will
// lose the track fields in the holes, so information from the holes has to be
// copied beforehand
   if (fNtracks == 0 || fCompact) return 0;
   fCompact = kTRUE;
   Int_t firsthole = fHoles.FirstSetBit();
   while (firsthole<fNtracks) {
      Int_t lastactive = fHoles.LastNullBit(fNtracks-1);
      if (lastactive < fNtracks) {
         // move last holes (if any)
         if (moveto && (fNtracks-lastactive-1>0)) moveto->AddTracks(*this, lastactive+1, fNtracks-1);
         fNtracks = lastactive+1;
         if (firsthole==fNtracks) return fNtracks;
      } else {
         // No active tracks left. First copy the hole track to the output
         if (moveto) moveto->AddTracks(*this, firsthole, firsthole+fNtracks-1);
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
   fSelected.ResetAllBits();
   fNselected = 0;
   return fNtracks;
}

//______________________________________________________________________________
Int_t GeantTrack_v::Reshuffle()
{
    //Printf("Reshuffling tracks on %p",this);
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
      fSelected.SetBitNumber(lastsel, false);
      firsthole = fSelected.FirstNullBit(firsthole+1);
      fNselected--;
   }
   return fNselected;
}

//______________________________________________________________________________
Bool_t GeantTrack_v::Contains(Int_t evstart, Int_t nevents) const
{
// Check if the array contains tracks from a given event range
   Int_t evend = evstart+nevents;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (fEventV[itr]>=evstart && fEventV[itr]<evend) return kTRUE;
   }
   return kFALSE;
}

//______________________________________________________________________________
void GeantTrack_v::Clear(Option_t *)
{
  //  Printf("Clearing tracks on %p",this);
    // Clear track content and selections
   fNselected = 0;
   fHoles.ResetAllBits();
   fSelected.ResetAllBits();
   fCompact = kTRUE;
   fNtracks = 0;
#ifdef __STAT_DEBUG_TRK
   fStat.Reset();
#endif

}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateStraight(Int_t ntracks, Double_t const *crtstep)
{
// Propagate first ntracks along a straight line (neutral particles, no mag.
// field or for last tiny step). The array has to be reshuffled after selecting
// the interesting particles using Select method.
// The crossing tracks get masked as holes in the array.

// Find next volume
   Printf("Propagate Straight %p", this);
   Int_t icrossed = 0;
   for (Int_t i=0; i<ntracks; i++) {
      if (fFrombdrV[i]) {
          Printf("Assigning next volume for track %d");
         *fPathV[i] = *fNextpathV[i];
         fStatusV[i] = kBoundary;
         icrossed++;
      }
   }

    for (Int_t i=0; i<ntracks; i++) {
      fPstepV[i] -= crtstep[i];
      assert( fPstepV[i] >= 0. && "STEP GOT NEGATIVE");
      fSafetyV[i] = 0;
      // Change path to reflect the physical volume for the current track; The
      // relevant path is fPath[i] if the frombdr flag is not set or fNextpath[i]
      // otherwise
      fXposV[i] += crtstep[i]*fXdirV[i];
      fYposV[i] += crtstep[i]*fYdirV[i];
      fZposV[i] += crtstep[i]*fZdirV[i];
      fNstepsV[i]++;
   }
   return icrossed;
}

//______________________________________________________________________________
void GeantTrack_v::PropagateInVolume(Int_t ntracks, const Double_t *crtstep)
{
#ifdef VERBOSE
	Printf("Propagate in Volume %d");
#endif
	// Propagate the selected tracks with crtstep values. The method is to be called
// only with  charged tracks in magnetic field. The tracks array has to be reshuffled.
   Double_t c = 0.;
   const Double_t *point = 0;
   const Double_t *newdir = 0;
   const Double_t bmag = gPropagator->fBmag;
   GeantThreadData *td = gPropagator->fThreadData[TGeoManager::ThreadId()];
   TGeoHelix *fieldp = td->fFieldPropagator;

   // TODO: to be vectorized
   for (Int_t i=0; i<ntracks; ++i) {
      // TODO: this function should just call/inline the next function here

	   // Reset relevant variables
      fFrombdrV[i] = kFALSE;
      fPstepV[i] -= crtstep[i];
      assert( fPstepV[i] >= 0. && "STEP GOT NEGATIVE");
      fSafetyV[i] -= crtstep[i];
      if (fSafetyV[i]<0.) fSafetyV[i] = 0.;
      fStepV[i] += crtstep[i];
      // Set curvature, charge
      c = std::fabs(kB2C*bmag/Pt(i));
// NOTE: vectorized treatment in TGeoHelix
      fieldp->SetXYcurvature(c);
      fieldp->SetCharge(fChargeV[i]);
      fieldp->SetHelixStep(std::fabs(TMath::TwoPi()*Pz(i)/(c*Pt(i))));
      fieldp->InitPoint(fXposV[i],fYposV[i],fZposV[i]);
      fieldp->InitDirection(fXdirV[i],fYdirV[i], fZdirV[i]);
      fieldp->UpdateHelix();
      fieldp->Step(crtstep[i]);
      point = fieldp->GetCurrentPoint();
      newdir = fieldp->GetCurrentDirection();
      fXposV[i] = point[0]; fYposV[i] = point[1]; fZposV[i] = point[2];
      fXdirV[i] = newdir[0]; fYdirV[i] = newdir[1]; fZdirV[i] = newdir[2];
   }
}

//______________________________________________________________________________
void GeantTrack_v::PropagateInVolumeSingle(Int_t i, Double_t crtstep)
{
#ifdef VERBOSE
	Printf("Propagate in Volume single %p track %d and step %lf", this, i, crtstep);
#endif
	// Propagate the selected tracks with crtstep values. The method is to be called
// only with  charged tracks in magnetic field. The tracks array has to be reshuffled.
   Double_t c = 0.;
   const Double_t *point = 0;
   const Double_t *newdir = 0;
   const Double_t bmag = gPropagator->fBmag;
   // Question: What has geometry to do with treading?
   GeantThreadData *td = gPropagator->fThreadData[TGeoManager::ThreadId()];

   TGeoHelix *fieldp = td->fFieldPropagator;
   // Reset relevant variables
   fFrombdrV[i] = kFALSE;
   fPstepV[i] -= crtstep;
#ifdef VERBOSE
   if( fPstepV[i] < 0 ) { Warning("PropagateinVolumeSingle", "PStep got negative %lf",fPstepV[i]);}
#endif
   //assert( fPstepV[i] >= 0. && "STEP GOT NEGATIVE");
   fSafetyV[i] -= crtstep;
   if (fSafetyV[i]<0.) fSafetyV[i] = 0.;
   fStepV[i] += crtstep;
   // Set curvature, charge
   c = std::fabs(kB2C*bmag/Pt(i));
// NOTE: vectorized treatment in TGeoHelix
   fieldp->SetXYcurvature(c);
   fieldp->SetCharge(fChargeV[i]);
   fieldp->SetHelixStep(std::fabs(TMath::TwoPi()*Pz(i)/(c*Pt(i))));
   fieldp->InitPoint(fXposV[i],fYposV[i],fZposV[i]);
   fieldp->InitDirection(fXdirV[i],fYdirV[i], fZdirV[i]);
   fieldp->UpdateHelix();
   fieldp->Step(crtstep);
   point = fieldp->GetCurrentPoint();
   newdir = fieldp->GetCurrentDirection();
   fXposV[i] = point[0]; fYposV[i] = point[1]; fZposV[i] = point[2];
   fXdirV[i] = newdir[0]; fYdirV[i] = newdir[1]; fZdirV[i] = newdir[2];
}

#ifdef USE_VECGEOM_NAVIGATOR

// provide here implementation calling VecGeom
void GeantTrack_v::NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep,
                       const Double_t *x, const Double_t *y, const Double_t *z,
                       const Double_t *dirx, const Double_t *diry, const Double_t *dirz,
                       VolumePath_t **pathin, VolumePath_t **pathout,
                       Double_t *step, Double_t *safe, Bool_t *isonbdr, const GeantTrack_v *trk)
{
   // Printf("In vec find next boundary and step\n");
    using vecgeom::SimpleNavigator;
    using vecgeom::Precision;
    using vecgeom::Vector3D;
    using vecgeom::GeoManager;
    typedef Vector3D<Precision> Vector3D_t;

    // VolumePath_t * a = new VolumePath_t( GeoManager::Instance().getMaxDepth() );

    SimpleNavigator nav;
    for (Int_t i=0; i<ntracks; ++i)
    {
#ifdef VERBOSE
      if( pstep[i] < 0. )
       {
        std::cerr << " NEGATIVE PSTEP " << pstep[i] << "\n";
        }
#endif

//    	a->Clear();
//    	nav.LocatePoint( GeoManager::Instance().world(),
//    			Vector3D_t( x[i], y[i], z[i] ), *a, true );
//        if( a->Top() != NULL && a->Top() != pathin[i]->Top() )
//         {
//             Printf("INCONSISTENT PATH, boundary state %d", isonbdr[i] );
//             a->GetCurrentNode()->Print();
//             pathin[i]->GetCurrentNode()->Print();
//             Printf("environment supposed path" );
//             nav.InspectEnvironmentForPointAndDirection(
//                             Vector3D_t( x[i], y[i], z[i] )  /*global pos */,
//                             Vector3D_t( dirx[i], diry[i], dirz[i] )  /*global dir*/ ,
//                             *pathin[i]);
//             Printf( "environment reported path" );
//             nav.InspectEnvironmentForPointAndDirection(
//                                         Vector3D_t( x[i], y[i], z[i] )  /*global pos*/ ,
//                                         Vector3D_t( dirx[i], diry[i], dirz[i] )  /*global dir*/ ,
//                                         *a);
//         }


	//        assert( a->Top() == pathin[i]->Top() );
    	nav.FindNextBoundaryAndStep( Vector3D_t( x[i], y[i], z[i] ) /* global pos */,
    			Vector3D_t( dirx[i], diry[i], dirz[i] ) /* global dir */,
    			*pathin[i], *pathout[i] /* the paths */,
    			TMath::Min(1.E20, pstep[i]),
    			step[i] );
        step[i] = TMath::Max(2.*gTolerance, step[i]);
        safe[i] = ( isonbdr[i] )? 0 : nav.GetSafety( Vector3D_t( x[i], y[i], z[i] ), *pathin[i] );
        safe[i] = (safe[i] < 0)? 0. : safe[i];

#ifdef CROSSCHECK
	//************
	// CROSS CHECK USING TGEO
	//************
	TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
	rootnav->ResetState();
	rootnav->SetCurrentPoint(x[i], y[i], z[i]);
	rootnav->SetCurrentDirection(dirx[i], diry[i], dirz[i]);
	TGeoBranchArray * tmp = pathin[i]->ToTGeoBranchArray();
	tmp->UpdateNavigator(rootnav);
	delete tmp;
	rootnav->FindNextBoundaryAndStep(TMath::Min(1.E20, pstep[i]), !isonbdr[i]);
	double stepcmp = TMath::Max(2*gTolerance,rootnav->GetStep());
	double safecmp = rootnav->GetSafeDistance();
	Printf("## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[i], step[i], stepcmp);
	Printf("## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf", pstep[i], isonbdr[i], safe[i], rootnav->Safety());

	// check nextpath
	tmp = pathout[i]->ToTGeoBranchArray();
	tmp->InitFromNavigator( rootnav );
	Printf("## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", pathout[i]->GetCurrentNode(), tmp->GetCurrentNode());
	Printf("## boundary result %d versus %d", pathout[i]->IsOnBoundary(), rootnav->IsOnBoundary());


	if( safe[i] != safecmp )
	{
		 nav.InspectSafetyForPoint(
		                                         Vector3D_t( x[i], y[i], z[i] )  /*global pos*/,
		                                         *pathin[i] );
		 //assert( ! "HELLO");
		 // correcting the safety to ROOT safety
		 step[i] = TMath::Max(2.*gTolerance, safecmp);
	}
#endif

	// onboundary with respect to new point
	isonbdr[i] = pathout[i]->IsOnBoundary();
#ifdef VERBOSE
	Printf("navfindbound on %p track %d with pstep %lf yields step %lf and safety %lf\n", this, i, pstep[i], step[i], safe[i]);
#endif
    }
 //   delete a;
}

#else
//______________________________________________________________________________
void GeantTrack_v::NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep,
                       const Double_t *x, const Double_t *y, const Double_t *z,
                       const Double_t *dirx, const Double_t *diry, const Double_t *dirz,
                       TGeoBranchArray **pathin, TGeoBranchArray **pathout,
                       Double_t *step, Double_t *safe, Bool_t *isonbdr, const GeantTrack_v *trk)
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
    //  assert( pstep[i]>= 0 );
	  nav->ResetState();
      nav->SetCurrentPoint(x[i], y[i], z[i]);
      nav->SetCurrentDirection(dirx[i], diry[i], dirz[i]);
 //     if( nav->FindNode( x[i], y[i], z[i] ) != pathin[i]->GetCurrentNode() )
 //     {
  //        Printf("INCONSISTENT PATH; boundarystatus %d, %s vs found %s", isonbdr[i],
   //               pathin[i]->GetCurrentNode()->GetName(),
    //              nav->FindNode( x[i], y[i], z[i] )->GetName());
    //  }
      pathin[i]->UpdateNavigator(nav);
      nav->FindNextBoundaryAndStep(TMath::Min(1.E20, pstep[i]), !isonbdr[i]);
      step[i] = TMath::Max(2*gTolerance,nav->GetStep());
      safe[i] = nav->GetSafeDistance();
#ifdef VERBOSE
      double bruteforces = nav->Safety();
      Printf("##TGEOM  ## TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg BRUTEFORCES %lg TOBOUND %d",
             i, isonbdr[i], pstep[i], step[i], safe[i], bruteforces, nav->IsOnBoundary());
#endif
      pathout[i]->InitFromNavigator(nav);

     // assert( safe[i]<=bruteforces );

#ifdef CROSSCHECK
      // crosscheck with what VECGEOM WOULD GIVE IN THIS SITUATION
      // ---------------------------------------------------------
      vecgeom::NavigationState vecgeom_in_state( vecgeom::GeoManager::Instance().getMaxDepth() );
      vecgeom::NavigationState vecgeom_out_state( vecgeom::GeoManager::Instance().getMaxDepth() );
      vecgeom_in_state = *pathin[i];
      vecgeom::SimpleNavigator vecnav;
      double vecgeom_step;
      typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;
      vecnav.FindNextBoundaryAndStep( Vector3D_t( x[i], y[i], z[i] ) /* global pos */,
                    Vector3D_t( dirx[i], diry[i], dirz[i] ) /* global dir */,
                    vecgeom_in_state, vecgeom_out_state /* the paths */,
                    TMath::Min(1.E20, pstep[i]),
                    vecgeom_step );
      vecgeom_step = TMath::Max(2*gTolerance, vecgeom_step);
      double vecgeom_safety;
      vecgeom_safety = vecnav.GetSafety( Vector3D_t( x[i], y[i], z[i] ),
                       vecgeom_in_state);
      vecgeom_safety = (vecgeom_safety < 0)? 0. : vecgeom_safety;
      Printf("--VECGEOM-- TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg TOBOUND %d",
                 i, isonbdr[i], pstep[i], vecgeom_step, vecgeom_safety, vecgeom_out_state.IsOnBoundary());
      // end crosscheck with what VECGEOM WOULD GIVE IN THIS SITUATION
      // ---------------------------------------------------------
#endif

      isonbdr[i] = nav->IsOnBoundary();
   }
}
#endif

#ifdef USE_VECGEOM_NAVIGATOR

/**
 * documentation for this function:
 * This function modifies the end path
 *
 */
void GeantTrack_v::NavIsSameLocation(Int_t ntracks, VolumePath_t ** start, VolumePath_t ** end, Bool_t *same)
{
#ifdef VERBOSE
    Printf("In NavIsSameLocation VECTOR MODE\n");
#endif
    // TODO: We should provide this function as a static function
    // TODO: use direction ( if needed )
    for(Int_t itr=0; itr<ntracks; ++itr)
    {
       NavIsSameLocationSingle(itr,start,end);
    }
}

#else
//______________________________________________________________________________
void GeantTrack_v::NavIsSameLocation(Int_t ntracks, VolumePath_t **start, VolumePath_t **end, Bool_t *same)
{
// Implementation of TGeoNavigator::IsSameLocation with vector input
// probably not used
   Printf("In NavIsSameLocation VECTOR MODE");
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   for (Int_t i=0; i<ntracks; i++) {
      nav->ResetState();
      nav->SetLastSafetyForPoint(0,0,0,0);
      nav->SetCurrentPoint(fXposV[i], fYposV[i], fZposV[i]);
      nav->SetCurrentDirection(fXdirV[i], fYdirV[i], fZdirV[i]);
      start[i]->UpdateNavigator(nav);
      same[i] = nav->IsSameLocation(fXposV[i],fYposV[i],fZposV[i],kTRUE);
      if (!same[i]) end[i]->InitFromNavigator(nav);
   }
}
#endif


#ifdef USE_VECGEOM_NAVIGATOR
// provide here implementation calling VecGeom
// what is this function doing?
// reseting navigation state flags ( about status of particle navigation: entering, leaving, etc... )
// setting last point to 0,0,0 and its safety to 0.?
// set current point and direction of track to navigator
// then update navigator caches with information from volumepath ( probably calculating the global matrix and updating local points in navigator ??)
//______________________________________________________________________________
Bool_t GeantTrack_v::NavIsSameLocationSingle(Int_t itr, VolumePath_t ** start,  VolumePath_t ** end)
{
#ifdef VERBOSE
     Printf("In NavIsSameLocation single %p for track %d", this, itr );
#endif
     // TODO: We should provide this function as a static function
    vecgeom::SimpleNavigator simplenav;
    vecgeom::NavigationState tmpstate( *end[itr] );

    // cross check with answer from ROOT
#ifdef CROSSCHECK
    TGeoBranchArray *sb = start[itr]->ToTGeoBranchArray();
    TGeoBranchArray *eb = end[itr]->ToTGeoBranchArray();
#endif

    // TODO: not using the direction yet here !!
    bool samepath = simplenav.HasSamePath(
                     vecgeom::Vector3D<vecgeom::Precision>(fXposV[itr],fYposV[itr],fZposV[itr]),
                     *start[itr],
                     tmpstate);
    if( !samepath )
    {
       *end[itr]=tmpstate;
    }

#ifdef CROSSCHECK
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0,0,0,0);
    nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
    nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
    sb->UpdateNavigator(nav);
    bool rootsame = nav->IsSameLocation(fXposV[itr],fYposV[itr],fZposV[itr],kTRUE);
    if( rootsame != samepath )
    {
        Printf("INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
        std::cout << vecgeom::Vector3D<vecgeom::Precision>(fXposV[itr],fYposV[itr],fZposV[itr]) << "\n";
        Printf("old state");sb->Print();
        nav->ResetState();
        nav->SetLastSafetyForPoint(0,0,0,0);
        nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
        nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
        sb->UpdateNavigator(nav);
        nav->InspectState();
        bool rootsame = nav->IsSameLocation(fXposV[itr],fYposV[itr],fZposV[itr],kTRUE);
        nav->InspectState();
        eb->InitFromNavigator(nav);
        Printf("new state");eb->Print();
        Printf("VERSUS VECGEOM OLD AND NEW");
        start[itr]->printVolumePath();
        end[itr]->printVolumePath();
    }
    else
    {
        Printf("CONSISTENT SAME LOCATION");
    }

    delete sb;
    delete eb;
#endif // CROSSCHECK
    return samepath;
}
#else
//______________________________________________________________________________
Bool_t GeantTrack_v::NavIsSameLocationSingle(Int_t itr, TGeoBranchArray **start, TGeoBranchArray **end)
{
// Implementation of TGeoNavigator::IsSameLocation for single particle
#ifdef VERBOSE
	Printf("IN Nav is same location SINGLE");
#endif
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0,0,0,0);
  nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
  nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  start[itr]->UpdateNavigator(nav);

#ifdef CROSSCHECK
  //* find out what VecGeom would do
   vecgeom::NavigationState vecgeom_in_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom::NavigationState vecgeom_out_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom_in_state = *start[itr];
   vecgeom::SimpleNavigator vecnav;
   typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;
   bool vecgeomresult = vecnav.HasSamePath( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ) /* global pos */,
                        vecgeom_in_state, vecgeom_out_state /* the paths */);
#endif
  if ( ! nav->IsSameLocation(fXposV[itr],fYposV[itr],fZposV[itr],kTRUE) ) {
    // this means that position is not in path given by path
    // the track has made a crossing --> recalculate where path went
   // most often this will be the same result as already present in end navigation state
#ifdef VERBOSE
     Printf("CORRECTING END STATE");
#endif
     end[itr]->InitFromNavigator(nav);
#ifdef CROSSCHECK
     if( vecgeomresult == true )
     {
        Warning("", "Problem with IsSameLocation");
        InspectIsSameLocation( itr );
     }
#endif
     return kFALSE;
  }
#ifdef CROSSCHECK
  if( vecgeomresult == false )
    {
       Warning("", "Problem with IsSameLocation");
       InspectIsSameLocation( itr );
    }
#endif
  return kTRUE;
}
#endif

#ifdef USE_VECGEOM_NAVIGATOR

//______________________________________________________________________________
void GeantTrack_v::PropagateBack(Int_t itr, Double_t crtstep)
{
#ifdef VERBOSE
    Printf("In propagateback %p for track %d\n", this, itr);
    fPathV[itr]->printVolumePath();
    fNextpathV[itr]->printVolumePath();
#endif
    typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;

    // This method is to be called after a successful crossing made by PropagateInField
    // to put the particle at the entry point.
    // check if particle is outside detector

    Bool_t outside = fNextpathV[itr]->IsOutside();


    // try to find out whether particle is entering where??
    Bool_t entering = kTRUE;
    // is the level of next path smaller than level of where it was coming from
    // this means that propagatin ba


//    Int_t level = nav->getCurrentLevel();
    Int_t level = fNextpathV[itr]->GetCurrentLevel();
    if ( level < fPathV[itr]->GetCurrentLevel() && !outside )
    {
         for ( Int_t lev=0; lev<=level; ++lev ){

            // if at some level higher up the path diverges do not do anything
            // this is something that should be implemented in the navigation state API
            if ( fNextpathV[itr]->At(lev) != fPathV[itr]->At(lev) ) break;
            if ( lev == level ) entering = kFALSE;
        }
    }


    // convert global coordinates to local coordinates
    Vector3D_t lpos, ldir;
    vecgeom::VPlacedVolume const * vol;
    if( ! outside )
    {
        vecgeom::TransformationMatrix transf = fNextpathV[itr]->TopMatrix();
        lpos = transf.Transform( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ));
        // transform of reversed direction ( because we are propagating back )
        ldir = transf.TransformDirection( Vector3D_t( -fXdirV[itr], -fYdirV[itr], -fZdirV[itr] ));

        // Put the navigation to the propagation location stored as next path
        vol = fNextpathV[itr]->Top();
    }
    else
    {
        lpos =  Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] );
        ldir =  Vector3D_t( -fXdirV[itr], -fYdirV[itr], -fZdirV[itr] );
       // Put the navigation to the propagation location stored as next path
        vol = vecgeom::GeoManager::Instance().world();
    }
    Double_t delta;
    if (entering) // go back to mother volume
       {
        assert( vol != NULL );
        // no matrix needed here because the world is placed with identity matrix
#ifdef VERBOSE
        if(outside) Printf("ENTERING; FROM OUTSIDE");
        if(!outside) Printf("ENTERING; INSIDE");
#endif
        delta = (outside)? vol->DistanceToIn( lpos, ldir ) : vol->DistanceToOut( lpos, ldir );
        }
    else // go back to daughter volume
    {
#ifdef VERBOSE
       Printf("NOT ENTERING; DISTANCETOIN");
#endif
       // TODO: check this !!!
       delta = fPathV[itr]->Top()->DistanceToIn( lpos, ldir );
    }
#ifdef VERBOSE
    Printf("Delta : %lf", delta);
#endif
    /*
   if (useDebug && (debugTrk<0 || itr==debugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }
*/
   if (delta>crtstep) {
/*
      if (useDebug && (debugTrk<0 || itr==debugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }
*/
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;

   // Propagate back to boundary and store position/direction
   fXposV[itr] -= delta*fXdirV[itr];
   fYposV[itr] -= delta*fYdirV[itr];
   fZposV[itr] -= delta*fZdirV[itr];
}

#else
//______________________________________________________________________________
void GeantTrack_v::PropagateBack(Int_t itr, Double_t crtstep)
{
#ifdef VERBOSE
	Printf("In propagateback for track %d\n", itr);
    fPathV[itr]->Print();
    fNextpathV[itr]->Print();
#endif
    // This method is to be called after a successful crossing made by PropagateInField
    // to put the particle at the volume entry point (the track curvature may
    // produce a boundary crossing macroscopically different than the expected one.
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    // Put the navigation to the propagation location stored as next path
    fNextpathV[itr]->UpdateNavigator(nav);
    // Backup the starting node location
    TGeoNode *checked = nav->GetCurrentNode();


   TGeoVolume *vol = checked->GetVolume();
   Double_t dir[3], ldir[3], ld[3];
   Double_t local[3], lp[3];
   Double_t delta;

   Bool_t outside = fNextpathV[itr]->IsOutside();

   // Swap track direction and compute distance back to boundary
   dir[0] = -fXdirV[itr]; dir[1] = -fYdirV[itr]; dir[2] = -fZdirV[itr];
   // Backup current depth in the hierarchy
   Int_t level = nav->GetLevel();
   Bool_t entering = kTRUE;
   // We need to check if:
   // 1. the track exited to a mother volume, in which case fNextpath should be
   // "contained" in fPath
   // - e.g. path = /a/b/c/d   nextpath = /a/b/c
   // 2. the track is entering a volume that did not containe the point before
   // crossing - in which case path and nextpath should diverge at a
   // level<current level
   // - e.g. path = /a/b/c     nextpath = /a/d
   if (level < fPathV[itr]->GetLevel() && !outside) {
      for (Int_t lev=0; lev<=level; lev++) {
         if (fNextpathV[itr]->GetNode(lev) != fPathV[itr]->GetNode(lev)) break;
         if (lev==level) entering = kFALSE;
      }
   }
   Double_t pos[3];
   pos[0] = fXposV[itr]; pos[1] = fYposV[itr]; pos[2] = fZposV[itr];
   // Convert to local frame
   fNextpathV[itr]->GetMatrix()->MasterToLocal(pos,local);
   fNextpathV[itr]->GetMatrix()->MasterToLocalVect(dir, ldir);
   // Depending if we enter or exit, we have to call the appropriate DistTo
   // method
   if (entering) {
      if (outside) delta = vol->GetShape()->DistFromOutside(local,ldir,3);
      else         delta = vol->GetShape()->DistFromInside(local,ldir,3);
   } else {
      checked = fPathV[itr]->GetNode(level+1);
      checked->MasterToLocal(local,lp);
      checked->MasterToLocalVect(ldir,ld);
      delta = checked->GetVolume()->GetShape()->DistFromOutside(lp,ld,3);
   }
#ifdef VERBOSE
   Printf("Delta : %lf", delta);
#endif

/*
   if (useDebug && (debugTrk<0 || itr==debugTrk)) {
      if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
      else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
   }
*/
   if (delta>crtstep) {
   // This signals sort of an error since the back propagation is larger than
  // the last step made for crossing the boundary
/*
      if (useDebug && (debugTrk<0 || itr==debugTrk)) {
         if (entering) Printf("   field-> track %d entering %s  at (%19.15f, %19.15f, %19.15f) crtstep=%19.15f delta=%19.15f", itr, vol->GetName(), xpos, ypos, zpos,crtstep, delta);
         else          Printf("   field-> track %d exiting %s, entering %s  crtstep=%19.15f delta=%19.15f", itr, checked->GetName(), vol->GetName(), crtstep, delta);
         Printf("Error propagating back %19.15f (forward %19.15f) track for track %d", delta, crtstep, itr);
         if (entering) Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", vol->GetName(), local[0],local[1],local[2],ldir[0],ldir[1],ldir[2]);
         else          Printf("%s Local: (%19.15f, %19.15f, %19.15f), ldir: (%19.15f, %19.15f, %19.15f)", checked->GetName(), lp[0],lp[1],lp[2],ld[0],ld[1],ld[2]);
      }
*/
      delta = 10*gTolerance;
   }
   delta -= 10*gTolerance;
   // Propagate back to boundary and store position/direction
   fXposV[itr] += delta*dir[0];
   fYposV[itr] += delta*dir[1];
   fZposV[itr] += delta*dir[2];
}
#endif

//______________________________________________________________________________
//NOTE(SW): this function is currently not called
Int_t GeantTrack_v::PropagateInField(Int_t ntracks, const Double_t *crtstep)
{
// Propagate with crtstep using the helix propagator. Mark crossing tracks as holes.
   Int_t icrossed = 0;
   Bool_t *same = new Bool_t[ntracks];
   PropagateInVolume(ntracks, crtstep);
   NavIsSameLocation(ntracks, fPathV, fNextpathV, same);
   for (Int_t itr=0; itr<ntracks; itr++) {
      if (same[itr]) continue;
      // Boundary crossed
      PropagateBack(itr, crtstep[itr]);
      // Update current path
      *fPathV[itr]=*fNextpathV[itr];
      fStatusV[itr] = kBoundary;
      fFrombdrV[itr] = kTRUE;
      MarkRemoved(itr);
      icrossed++;
   }
   return icrossed;
}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateInFieldSingle(Int_t itr, Double_t crtstep, Bool_t checkcross)
{
// Propagate with crtstep using the helix propagator. Mark crossing tracks as holes.
#ifdef VERBOSE
	Printf("propagate in field Single with crtstep %lf and checkcross %d",crtstep, checkcross);
#endif
   PropagateInVolumeSingle(itr, crtstep);
#ifdef VERBOSE
   bool isnullstep = crtstep < 1E-8;
   if( isnullstep )
   {
      Warning("PROPAGATEINFIELDSINGLE"," STEP IS NULL ");
      if(! fPathV[itr]->IsOutside() ) fPathV[itr]->GetCurrentNode()->Print();
      if(! fNextpathV[itr]->IsOutside() ) fNextpathV[itr]->GetCurrentNode()->Print();
   }
#endif
   if (checkcross && ! NavIsSameLocationSingle(itr, fPathV, fNextpathV)) {
#ifdef VERBOSE
     if(isnullstep){
     Warning("PROPAGATEINFIELDSINGLE","WOW DESPITE NULL STEP WE ARE HERE !!!");
      }
      // Boundary crossed
      if( fPathV[itr]->GetCurrentNode() == fNextpathV[itr]->GetCurrentNode() )
      {
         Warning("PropagateInFieldSingle","Crossed boundary but actually the paths are the same %d", itr);
         fPathV[itr]->Print();
         fNextpathV[itr]->Print();
         Printf("%lf, %lf, %lf\n", fXposV[itr], fYposV[itr], fZposV[itr]);
         // let me see what the navigator gives me:
         TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
         nav->FindNode( fXposV[itr], fYposV[itr], fZposV[itr] );
         nav->InspectState();
         Warning("PropagateInFieldSingle","END OF WARNING MESSAGE %d", itr);
      }
#endif
     PropagateBack(itr, crtstep);

#ifdef CROSSCHECK
     // just to be sure: check again consistency
     TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
     nav->FindNode( fXposV[itr], fYposV[itr], fZposV[itr] );
     if( nav->GetCurrentNode() != fNextpathV[itr]->GetCurrentNode() )
     {
       Warning("PropagateInFieldSingle","PropagateBackWasNotConsistent");
       nav->GetCurrentNode()->Print();
#ifdef USE_VECGEOM_NAVIGATOR
       using vecgeom::SimpleNavigator;
       using vecgeom::NavigationState;
       using vecgeom::Precision;
       using vecgeom::Vector3D;
       using vecgeom::GeoManager;
       typedef Vector3D<Precision> Vector3D_t;
       SimpleNavigator nav;
       NavigationState state( GeoManager::Instance().getMaxDepth() );
       nav.LocatePoint( GeoManager::Instance().world(),
       Vector3D<Precision>(fXposV[itr], fYposV[itr], fZposV[itr]),
       state, true );
       // state.GetCurrentNode()->Print();
#endif
       Warning("PropagateInFieldSingle","END OF WARNING PropagateBackWasNotConsistent");
     }
#endif

     // Update current path
     *fPathV[itr]=*fNextpathV[itr];
     fStatusV[itr] = kBoundary;
     fFrombdrV[itr] = kTRUE;
     MarkRemoved(itr);
     return 1;
   }
#ifdef VERBOSE
   else
   {
     if( isnullstep && checkcross )
      {
       Warning("PropagateInFieldSingle","NULLSTEP SHOULD CROSS BUT DID NOT with COUNTER STATUS %d", fIzeroV[itr]);
      }
   }
#endif
   return 0;
}

//______________________________________________________________________________
Int_t GeantTrack_v::SortByStatus(TrackStatus_t status)
{
// Sort tracks by a given status.
   Int_t nsel = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (fStatusV[itr] == status) {
         Select(itr);
         nsel++;
      }
   }
   if (nsel) Reshuffle();
   return nsel;
}

//______________________________________________________________________________
Int_t GeantTrack_v::RemoveByStatus(TrackStatus_t status, GeantTrack_v &output)
{
// Remove tracks with given status from the container to the output vector,
// then compact.
   Int_t nremoved = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (fStatusV[itr] == status) {
         MarkRemoved(itr);
         nremoved++;
      }
   }
   if (!fCompact) Compact(&output);
   return nremoved;
}

//______________________________________________________________________________
void GeantTrack_v::PrintTrack(Int_t itr) const
{
// Print info for a given track
      const char* status[7] = {"alive", "killed", "boundary", "exitSetup", "physics","postponed","new"};
#ifdef USE_VECGEOM_NAVIGATOR
      printf("Object %p, Track %d: evt=%d slt=%d part=%d pdg=%d g5c=%d chg=%d proc=%d izr=%d nstp=%d spc=%d status=%s mass=%g\
              xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g pstp=%g stp=%g snxt=%g saf=%g bdr=%d\n\n",
              this, itr, fEventV[itr],fEvslotV[itr], fParticleV[itr], fPDGV[itr],
              fG5codeV[itr], fChargeV[itr], fProcessV[itr],fIzeroV[itr],fNstepsV[itr],
              (Int_t)fSpeciesV[itr], status[Int_t(fStatusV[itr])], fMassV[itr], fXposV[itr],fYposV[itr],fZposV[itr],
              fXdirV[itr],fYdirV[itr],fZdirV[itr],fPV[itr],fEV[itr],fPstepV[itr], fStepV[itr], fSnextV[itr],fSafetyV[itr],
              fFrombdrV[itr]);

#else
      TString path;
      fPathV[itr]->GetPath(path);
      TString nextpath;
      fNextpathV[itr]->GetPath(nextpath);
      printf("Track %d: evt=%d slt=%d part=%d pdg=%d g5c=%d chg=%d proc=%d izr=%d nstp=%d spc=%d status=%s mass=%g xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g pstp=%g stp=%g snxt=%g saf=%g bdr=%d\n pth=%s npth=%s\n",
              itr, fEventV[itr],fEvslotV[itr], fParticleV[itr], fPDGV[itr],
              fG5codeV[itr], fChargeV[itr], fProcessV[itr],fIzeroV[itr],fNstepsV[itr],
              (Int_t)fSpeciesV[itr], status[Int_t(fStatusV[itr])], fMassV[itr], fXposV[itr],fYposV[itr],fZposV[itr],
              fXdirV[itr],fYdirV[itr],fZdirV[itr],fPV[itr],fEV[itr],fPstepV[itr], fStepV[itr], fSnextV[itr],fSafetyV[itr],
              fFrombdrV[itr], path.Data(), nextpath.Data());
#endif
}


//______________________________________________________________________________
void GeantTrack_v::PrintTracks() const
{
// Print all tracks
   for (Int_t i=0; i<fNtracks; i++) PrintTrack(i);
}

#ifdef USE_VECGEOM_NAVIGATOR
void GeantTrack_v::ComputeTransportLength(Int_t ntracks)
{
#ifdef VERBOSE
	Printf("Compute Transport Length");
#endif
	static Int_t icalls = 0;
    icalls++;
    Int_t itr;
    // call the vector interface of GeantTrack_v
    NavFindNextBoundaryAndStep(ntracks, fPstepV, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV,
                                 fPathV, fNextpathV, fSnextV, fSafetyV, fFrombdrV, this);
    // now we should have updated everything

    // perform a couple of additional checks/ set status flags and so on
    for (itr=0; itr<ntracks; ++itr) {
         if ((fNextpathV[itr]->IsOutside() && fSnextV[itr]<1.E-6) || fSnextV[itr]>1.E19) fStatusV[itr] = kExitingSetup;
         if (fFrombdrV[itr] && fSnextV[itr]<2.*gTolerance) {
            // Make sure track crossed
            fIzeroV[itr]++;
            if (fIzeroV[itr] > 10) {
               fStatusV[itr] = kKilled;
               Printf("### track %d had to be killed due to crossing problems", fParticleV[itr]);
               continue;
            }
            // assert(! "PLEASE IMPLEMENT");
            Warning("ComputeTransportLength","code tried to call TGEO navigator with potentially corrupted state");
            //nav->FindNextBoundaryAndStep(1.E30, kFALSE);
            //fNextpathV[itr]->InitFromNavigator(nav);
            //fSnextV[itr] += nav->GetStep();
         }
   //      if (fSnextV[itr]>2.*gTolerance) fIzeroV[itr] = 0;
      }
}
#else
//______________________________________________________________________________
void GeantTrack_v::ComputeTransportLength(Int_t ntracks)
{
#ifdef VERBOSE
	Printf("Compute Transport Length %p", this);
#endif
	// Computes snext and safety for an array of tracks. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.
   static Int_t icalls = 0;
   icalls++;
  for (Int_t itr=0; itr<ntracks; ++itr) {
        ComputeTransportLengthSingle(itr);
   }
   //---------------- alternative code ---------------------------------------
//   Int_t itr;
//   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
//   NavFindNextBoundaryAndStep(ntracks, fPstepV, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV,
//                              fPathV, fNextpathV, fSnextV, fSafetyV, fFrombdrV, this);
//   for (itr=0; itr<ntracks; itr++) {
//      if ((fNextpathV[itr]->IsOutside() && fSnextV[itr]<1.E-6) || fSnextV[itr]>1.E19) fStatusV[itr] = kExitingSetup;
//      if (fFrombdrV[itr] && fSnextV[itr]<=2.*gTolerance) {
//         // Make sure track crossed
//         fIzeroV[itr]++;
//         if (fIzeroV[itr] > 10) {
//            fStatusV[itr] = kKilled;
//            Printf("### track %d had to be killed due to crossing problems", fParticleV[itr]);
//            continue;
//         }
//         // what is the state of the navigator here?? very likely not the state of track itr
//         Warning("ComputeTransportLength","code tried to call TGEO navigator with potentially corrupted state");
//         nav->FindNextBoundaryAndStep(1.E30, kFALSE);
//         fNextpathV[itr]->InitFromNavigator(nav);
//         fSnextV[itr] += nav->GetStep();
//      }
////      if (fSnextV[itr]>2.*gTolerance) fIzeroV[itr] = 0;
//   }
}
#endif


#ifdef USE_VECGEOM_NAVIGATOR

void GeantTrack_v::ComputeTransportLengthSingle(Int_t itr)
{
// Computes snext and safety for a single track. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.
   static Int_t icalls = 0;
   icalls++;


      // inits navigator with current state
   //TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   //nav->ResetState();
   //nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
   //nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
   //fPathV[itr]->UpdateNavigator(nav);
   //nav->SetLastSafetyForPoint(fSafetyV[itr], fXposV[itr], fYposV[itr], fZposV[itr]);
   //nav->FindNextBoundaryAndStep( TMath::Min(1.E20, fPstepV[itr]), !fFrombdrV[itr] );

   //
   using vecgeom::SimpleNavigator;
   using vecgeom::Precision;
   using vecgeom::Vector3D;
   typedef Vector3D<Precision> Vector3D_t;

   vecgeom::SimpleNavigator nav;
   double step;
   nav.FindNextBoundaryAndStep(Vector3D_t(fXposV[itr], fYposV[itr], fZposV[itr]),
                               Vector3D_t(fXdirV[itr], fYdirV[itr], fZdirV[itr]),
                               *fPathV[itr],
                               *fNextpathV[itr],
                               TMath::Min(1.E20, fPstepV[itr]),
                               step);

   // get back step, safety, new geometry path, and other navigation information
   fSnextV[itr] = TMath::Max(2*gTolerance,step);
   fSafetyV[itr] = (fFrombdrV[itr])? 0 : nav.GetSafety( Vector3D_t(fXposV[itr],fYposV[itr],fZposV[itr]),
                                                        *fPathV[itr] );
   fSafetyV[itr] = (fSafetyV[itr]<0)? 0. : fSafetyV[itr];
   fFrombdrV[itr] = fNextpathV[itr]->IsOnBoundary();

   // if outside detector or enormous step mark particle as exiting the detector
   if (fNextpathV[itr]->IsOutside() || fSnextV[itr]>1.E19) fStatusV[itr] = kExitingSetup;

   // force track to cross under certain conditions
   if (fFrombdrV[itr] && fSnextV[itr] < 2.*gTolerance) {
      // Make sure track crossed
      fIzeroV[itr]++;
      if (fIzeroV[itr] > 10) {
         fStatusV[itr] = kKilled;
         Printf("### track %d had to be killed due to crossing problems", fParticleV[itr]);
         return;
      }
      // we have to make this general like above
      // this is moving the particle to the next volume??
      // nav.FindNextBoundaryAndStep(1.E30, kFALSE);
      // to be changed !!
      // fNextpathV[itr]->InitFromNavigator(nav);
      // to be changed !!
      // fSnextV[itr] += nav->GetStep();
      // assert( ! "NOT IMPLEMENTED" );
      printf("Funny situation; Finish implementation of ComputeTransportLength; ask Andrei!\n");
   }
   //   if (fSnextV[itr]>2.*gTolerance) fIzeroV[itr] = 0;
}

#else

void GeantTrack_v::ComputeTransportLengthSingle(Int_t itr)
{
// Computes snext and safety for a single track. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.
   static Int_t icalls = 0;
   static Int_t warnings=0;
   icalls++;
#ifdef VERBOSE
   Printf("Compute Transport Length Single %p on track", this, itr);
#endif
   // inits navigator with current state
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   nav->ResetState();
   nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
   nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
   fPathV[itr]->UpdateNavigator(nav);
  // nav->SetLastSafetyForPoint(fSafetyV[itr], fXposV[itr], fYposV[itr], fZposV[itr]);
   fSafetyV[itr] = nav->Safety();
   nav->FindNextBoundaryAndStep( TMath::Min(1.E20, fPstepV[itr]), !fFrombdrV[itr] );

   // get back step, safety, new geometry path, and other navigation information
   fSnextV[itr] = TMath::Max(2*gTolerance,nav->GetStep());

   fNextpathV[itr]->InitFromNavigator(nav);

#ifdef CROSSCHECK
   vecgeom::NavigationState vecgeom_in_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom::NavigationState vecgeom_out_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom_in_state = *fPathV[itr];
   vecgeom::SimpleNavigator vecnav;
   double vecgeom_step;
   typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;
   vecnav.FindNextBoundaryAndStep( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ) /* global pos */,
                       Vector3D_t( fXdirV[itr], fYdirV[itr], fZdirV[itr] ) /* global dir */,
                       vecgeom_in_state, vecgeom_out_state /* the paths */,
                       TMath::Min(1.E20, fPstepV[itr]),
                       vecgeom_step );
   vecgeom_step = TMath::Max(2*gTolerance, vecgeom_step);
   double vecgeom_safety;
   vecgeom_safety = vecnav.GetSafety( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ),
   vecgeom_in_state);
   vecgeom_safety = (vecgeom_safety < 0)? 0. : vecgeom_safety;
//   Printf("--VECGEOM-- TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg TOBOUND %d",
  //                  i, isonbdr[i], pstep[i], vecgeom_step, vecgeom_safety, vecgeom_out_state.IsOnBoundary());
   if( ! ( vecgeom_out_state.IsOutside() || fNextpathV[itr]->IsOutside() ) &&
        ( vecgeom_out_state.GetCurrentNode() != fNextpathV[itr]->GetCurrentNode() ) )
   {
	   warnings++;
	   InspectGeometryState(itr);
	   Warning("","NavigationWarning");
	   if(warnings>100) assert(! "NOT SAME RESULT FOR NEXTNODE" );
   }
#endif
   fFrombdrV[itr] = nav->IsOnBoundary();

   // if outside detector or enormous step mark particle as exiting the detector
   if (fNextpathV[itr]->IsOutside() || fSnextV[itr]>1.E19) fStatusV[itr] = kExitingSetup;

   // force track to cross under certain conditions
   if (fFrombdrV[itr] && fSnextV[itr]<2.*gTolerance) {
      // Make sure track crossed
      fIzeroV[itr]++;
      if (fIzeroV[itr] > 10) {
         fStatusV[itr] = kKilled;
         Printf("### track %d had to be killed due to crossing problems", fParticleV[itr]);
         return;
      }
      // try go give particle a push
      nav->FindNextBoundaryAndStep(1.E30, kFALSE);
      fNextpathV[itr]->InitFromNavigator(nav);
      fSnextV[itr] += nav->GetStep();
   }
//   if (fSnextV[itr]>2.*gTolerance) fIzeroV[itr] = 0;
}

#endif

void GeantTrack_v::InspectIsSameLocation( Int_t itr ) const
{
#ifndef USE_VECGEOM_NAVIGATOR
    Printf("STARTINSPECT----------------------------------");
    TGeoNavigator * nav = gGeoManager->GetCurrentNavigator();
    nav->ResetState();
   // nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
   // nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
    fPathV[itr]->UpdateNavigator(nav);

    bool rootissame = nav->IsSameLocation(fXposV[itr],fYposV[itr],fZposV[itr],kTRUE);

      //* find out what VecGeom would do
    vecgeom::NavigationState vecgeom_in_state( vecgeom::GeoManager::Instance().getMaxDepth() );
    vecgeom::NavigationState vecgeom_out_state( vecgeom::GeoManager::Instance().getMaxDepth() );
    vecgeom_in_state = *fPathV[itr];
    vecgeom::SimpleNavigator vecnav;
    typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;
    bool vecgeomresult = vecnav.HasSamePath( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ) /* global pos */,
                            vecgeom_in_state, vecgeom_out_state /* the paths */);

    fPathV[itr]->Print();
    vecgeom_in_state.printVolumePath();
    Printf("POINT( %lg, %lg, %lg), DIR( %lf, %lf, %lf)",
            fXposV[itr], fYposV[itr], fZposV[itr],
            fXdirV[itr], fYdirV[itr], fZdirV[itr]);
    Printf("ROOT: %d VecGeom %d", rootissame, vecgeomresult);
    TGeoBranchArray tmp( *fPathV[itr]);
    tmp.InitFromNavigator(nav);
    tmp.Print();
    vecgeom_out_state.printVolumePath();
    Printf("ENDINSPECT----------------------------------");
#endif
}

void GeantTrack_v::InspectGeometryState( Int_t itr ) const
{
#ifndef USE_VECGEOM_NAVIGATOR
// inits navigator with current state
   double rootsafety;
   double rootsnext;
   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
   TGeoBranchArray rootnextpath( *fPathV[itr] );
   nav->ResetState();
   nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
   nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
   fPathV[itr]->UpdateNavigator(nav);
  // nav->SetLastSafetyForPoint(fSafetyV[itr], fXposV[itr], fYposV[itr], fZposV[itr]);
   rootsafety = nav->Safety();
   nav->FindNextBoundaryAndStep( TMath::Min(1.E20, fPstepV[itr]), !fFrombdrV[itr] );
   // get back step, safety, new geometry path, and other navigation information
   rootsnext = TMath::Max(2*gTolerance,nav->GetStep());
   rootnextpath.InitFromNavigator(nav);

   //
   vecgeom::NavigationState vecgeom_in_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom::NavigationState vecgeom_out_state( vecgeom::GeoManager::Instance().getMaxDepth() );
   vecgeom_in_state = *fPathV[itr];
   vecgeom::SimpleNavigator vecnav;
   double vecgeom_step;
   typedef vecgeom::Vector3D<vecgeom::Precision> Vector3D_t;
   vecnav.FindNextBoundaryAndStep( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ) /* global pos */,
                       Vector3D_t( fXdirV[itr], fYdirV[itr], fZdirV[itr] ) /* global dir */,
                       vecgeom_in_state, vecgeom_out_state /* the paths */,
                       TMath::Min(1.E20, fPstepV[itr]),
                       vecgeom_step );
   vecgeom_step = TMath::Max(2*gTolerance, vecgeom_step);
   double vecgeom_safety;
   vecgeom_safety = vecnav.GetSafety( Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ),
   vecgeom_in_state);
   vecgeom_safety = (vecgeom_safety < 0)? 0. : vecgeom_safety;
fPathV[itr]->Print();
vecgeom_in_state.printVolumePath();
rootnextpath.Print();
Printf("ROOT: CURRENTNODEPOINTER %p", fPathV[itr]->GetCurrentNode());
Printf("ROOT: POINT( %lg, %lg, %lg), DIR( %lf, %lf, %lf), NORMPATH ( %lg )",
        fXposV[itr], fYposV[itr], fZposV[itr],
        fXdirV[itr], fYdirV[itr], fZdirV[itr],
        sqrt(fXdirV[itr]*fXdirV[itr] + fYdirV[itr]*fYdirV[itr] + fZdirV[itr]*fZdirV[itr])-1);
if( rootnextpath.GetCurrentNode()!=NULL )
{
Printf("ROOT: PSTEP %lf, SNEXT %lf, SAFETY %lf, NEXTNODE %p (level %d)",fPstepV[itr],rootsnext,rootsafety,
    rootnextpath.GetCurrentNode(), rootnextpath.GetLevel());
}
else
{
    Printf("ROOT: PSTEP %lf, SNEXT %lf, SAFETY %lf, NEXTNODE %p",fPstepV[itr],rootsnext,rootsafety,
    rootnextpath.GetCurrentNode());
}
if( vecgeom_out_state.GetCurrentNode()!=NULL )
{
Printf("VGEO: PSTEP %lf, SNEXT %lf, SAFETY %lf, NEXTNODE %p",fPstepV[itr],vecgeom_step,vecgeom_safety,
    vecgeom_out_state.GetCurrentNode());
}
else
{
Printf("VGEO: PSTEP %lf, SNEXT %lf, SAFETY %lf, NEXTNODE %p",fPstepV[itr],vecgeom_step,vecgeom_safety,
    vecgeom_out_state.GetCurrentNode());
}
// check if paths are actually consistent
  vecgeom::NavigationState vecgeom_check_state( vecgeom::GeoManager::Instance().getMaxDepth() );
  vecnav.LocatePoint( vecgeom::GeoManager::Instance().world(), Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ),
      vecgeom_check_state, true );
  Printf("VecGeom actually is here");
  vecgeom_check_state.printVolumePath();

  // if( vecgeom_in_state.Top() && vecgeom_check_state.Top() == vecgeom_in_state.Top() )
 // {
 //  vecgeom_check_state.printVolumePath();
 //  Printf("VecGeom PATH consistent");
 // }
 // else
 //{
 //  Printf("VecGEOM PATH WAS NOT CONSISTENT");
 // }

  // check if paths are actually consistent
  TGeoNode *tmp = nav->FindNode( fXposV[itr], fYposV[itr], fZposV[itr] );
  rootnextpath.InitFromNavigator(nav);
  Printf("ROOT actually is here");
  rootnextpath.Print();

  vecgeom::GeoManager::Instance().gVERBOSE=1;
  vecnav.InspectEnvironmentForPointAndDirection(
      Vector3D_t( fXposV[itr], fYposV[itr], fZposV[itr] ) /* global pos */,
                  Vector3D_t( fXdirV[itr], fYdirV[itr], fZdirV[itr] ),
      vecgeom_in_state );
  vecgeom::GeoManager::Instance().gVERBOSE=0;
#endif
}

//______________________________________________________________________________
TransportAction_t GeantTrack_v::PostponedAction() const
{
// Check the action to be taken according the current policy
   if (!fNtracks) return kDone;
   // Temporary hook
   if (fNtracks==1) {
//      if (gPropagator->GetPolicy()<GeantPropagator::kPrioritize)
//           return kPostpone;
//      else
           return kSingle;
   }
   return kVector;
}

//
/**
 * propagate tracks and return number of crossed tracks
 * builds an output track at same time
 */
//______________________________________________________________________________
Int_t GeantTrack_v::PropagateTracks(GeantTrack_v &output)
{
#ifdef VERBOSE
    Printf("In propagate tracks\n");
#endif
    // Propagate the ntracks in the current volume with their physics steps (already
    // computed)
    // Vectors are pushed downstream when efficient.

   // check consistency of this vector track: all particles should be in same logical volume
   bool insame=true;
   for( Int_t i=0;i<fNtracks-1; ++i )
   {
    if(fPathV[i]->GetCurrentNode()->GetVolume()
            != fPathV[i+1]->GetCurrentNode()->GetVolume()){
        Warning("","basket condition broken %s %s", fPathV[i]->GetCurrentNode()->GetVolume(),
                fPathV[i+1]->GetCurrentNode()->GetVolume());
        insame = false;
    }
   }
   assert( insame && "BASKET CONDITION BROKEN ");

   // Check if tracking the remaining tracks can be postponed
   TransportAction_t action = PostponedAction();
   if (action == kPostpone) {
      PostponeTracks(output);
      return 0;
   }
   if (action != kVector) return PropagateTracksSingle(output,0);

   // Compute transport length in geometry, limited by the physics step
   ComputeTransportLength(fNtracks);

   Int_t itr = 0;
   Int_t icrossed = 0;
   Int_t nsel = 0;
   Double_t c;
   const Double_t bmag = gPropagator->fBmag;

   // Remove dead tracks, propagate neutrals
   // TODO: should be vectorized
   for (itr=0; itr<fNtracks; itr++) {
      // Mark dead tracks for copy/removal
      if (fStatusV[itr] == kKilled) {
         MarkRemoved(itr);
         continue;
      }

      // Propagate straight tracks to the precomputed location and update state,
      // then mark them for copy/removal
      // (Inlined from PropagateStraight)
      if (fChargeV[itr]==0 || bmag<1.E-10) {
        // Do straight propagation to physics process or boundary
         Printf("Straight line propagation");
        if (fFrombdrV[itr]) {
            //SetPath(itr, *fNextpathV[itr]);
            *fPathV[itr] = *fNextpathV[itr];
            if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
            else                          fStatusV[itr] = kBoundary;
            icrossed++;
       } else {
            fStatusV[itr] = kPhysics;
       }
       fPstepV[itr] -= fSnextV[itr];
         assert( fPstepV[itr] >= 0. && "STEP GOT NEGATIVE");
         fStepV[itr] += fSnextV[itr];
         fSafetyV[itr] -= fSnextV[itr];
         if (fSafetyV[itr]<gTolerance) fSafetyV[itr] = 0;
         fXposV[itr] += fSnextV[itr]*fXdirV[itr];
         fYposV[itr] += fSnextV[itr]*fYdirV[itr];
         fZposV[itr] += fSnextV[itr]*fZdirV[itr];
         fNstepsV[itr]++;
         MarkRemoved(itr);
      }
   } // end first loop over tracks ( to treat neutrals )
   // Compact remaining tracks and move the removed oned to the output container
   if (!fCompact) Compact(&output);
   // Check if tracking the remaining tracks can be postponed
   action = PostponedAction();
   switch (action) {
      case kDone:
         return icrossed;
      case kSingle:
         icrossed += PropagateTracksSingle(output,1);
         return icrossed;
      case kPostpone:
         PostponeTracks(output);
         return icrossed;
      case kVector:
         break;
   }
   // Continue with vectorized mode ...

   // REMAINING ONLY CHARGED TRACKS IN MAG. FIELD
   Double_t *steps = new Double_t[fNtracks];
   // Select tracks that undergo safe steps to physics process
   nsel = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      if (fPstepV[itr]<fSafetyV[itr]) {
         Select(itr);
         fStatusV[itr] = kPhysics;
         nsel++;
      }
   }
   // fPstep array is contiguous after the last Reshuffle operation
   if (nsel) {
      Reshuffle();
      memcpy(steps, fPstepV, nsel*sizeof(Double_t));
      // This category of tracks can make the full step to the physics process
      PropagateInVolume(nsel, steps);
      // Move these tracks to the output container
      output.AddTracks(*this, 0, nsel-1);
      RemoveTracks(0, nsel-1);
   }
   delete [] steps;

   // CODE repetition ( same as above )
   action = PostponedAction();
   switch (action) {
      case kDone:
         return icrossed;
      case kSingle:
         icrossed += PropagateTracksSingle(output,1);
         return icrossed;
      case kPostpone:
         PostponeTracks(output);
         return icrossed;
      case kVector:
         break;
   }

   // Continue with vectorized mode ...
   // Check if we can propagate to boundary
   // those are the track that will hit a boundary

   nsel = 0;
   for (Int_t itr=0; itr<fNtracks; itr++) {
      c = Curvature(itr);
      if (0.25*c*fSnextV[itr]<1E-6 && fSnextV[itr]<1E-3 && fSnextV[itr]<fPstepV[itr]-1E-6) {
         // Propagate with snext and check if we crossed
         if (fIzeroV[itr]>10) fSnextV[itr] = 1.E-3;
         icrossed += PropagateInFieldSingle(itr, fSnextV[itr]+10*gTolerance, kTRUE);
         if (fSnextV[itr]<1.E-6) fIzeroV[itr]++;
         else fIzeroV[itr] = 0;
         // Crossing tracks have the correct status and are now marked for removal
         gPropagator->fNsnextSteps++;
         if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
         continue; // -> to next track
      }
      // Track has safety<pstep but next boundary not close enough.
      // We propagate in field with the safety value.
      // Protection in case we have a very small safety after entering a new volume
      const Double_t fs = 0.1;
//      if (fSafetyV[itr] < 1.E-6) {
//         fSafetyV[itr] = 1.E-3;
      if (/*fSafetyV[itr] < 0.001 &&*/ fSafetyV[itr] < fs*fSnextV[itr]) {
         // Track getting away from boundary. Work to be done here
         // In principle we need a safety value for the content of the current volume only
         // This does not behave well on corners...
         // ... so we peek a small value and chech if this crosses, than recompute safety
//
         Double_t delta = fs*fSnextV[itr];
         if (c*fSnextV[itr]<0.1) delta = 0.5*fSnextV[itr]*fSnextV[itr]*c;
         else if (c*fSnextV[itr]>10.) delta = TMath::Max(1/c, delta);
         fSafetyV[itr] = TMath::Max(fSafetyV[itr], delta);
         fIzeroV[itr]++;
         if (fIzeroV[itr] > 10) fSafetyV[itr] = 0.5*fSnextV[itr];
         icrossed += PropagateInFieldSingle(itr, fSafetyV[itr], kTRUE);
      } else {
         if (fIzeroV[itr] > 10) {
            // Propagate with snext
            icrossed += PropagateInFieldSingle(itr, fSnextV[itr]+10*gTolerance, kTRUE);
            fIzeroV[itr] = 0;
            gPropagator->fNsnextSteps++;
            if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
            continue; // -> to next track
         }
         if (fSafetyV[itr]<1.E-3) fIzeroV[itr]++;
         // Propagate with safety without checking crossing
         icrossed += PropagateInFieldSingle(itr, fSafetyV[itr], kFALSE);
      }
      gPropagator->fNsafeSteps++;
//      if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
   }
   // Compact remaining tracks and move the removed ones to the output container
   if (!fCompact) Compact(&output);
   // Remaining tracks have been partially propagated, they need to be
   // transported again after applying continuous energy loss
//   if (fNtracks) ComputeTransportLength(fNtracks);
   return icrossed;
}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateTracksSingle(GeantTrack_v &output, Int_t stage)
{
// Propagate the tracks with their selected steps in a single loop,
// starting from a given stage.
   Int_t itr = 0;
   Int_t icrossed = 0;
   Double_t step, c;
   const Double_t bmag = gPropagator->fBmag;
   for (itr=0; itr<fNtracks; itr++) {
      if (fStatusV[itr] == kKilled) {
         MarkRemoved(itr);
         continue;
      }
      // Compute transport length in geometry, limited by the physics step
      ComputeTransportLengthSingle(itr);
      // Stage 0: straight propagation
      if (stage==0) {
      // nice code repetition
         if (fChargeV[itr]==0 || bmag<1.E-10) {
         // Do straight propagation to physics process or boundary
             Printf("Straight line propagation (Single)");
             if (fFrombdrV[itr]) {
               *fPathV[itr] = *fNextpathV[itr];
               if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
               else                          fStatusV[itr] = kBoundary;
               icrossed++;
            } else {
               fStatusV[itr] = kPhysics;
            }
            fPstepV[itr] -= fSnextV[itr];
            assert( fPstepV[itr] >= 0. && "STEP GOT NEGATIVE");
            fSafetyV[itr] -= fSnextV[itr];
            if (fSafetyV[itr]<gTolerance) fSafetyV[itr] = 0;
            fXposV[itr] += fSnextV[itr]*fXdirV[itr];
            fYposV[itr] += fSnextV[itr]*fYdirV[itr];
            fZposV[itr] += fSnextV[itr]*fZdirV[itr];
            fNstepsV[itr]++;
            MarkRemoved(itr);
            continue;
         }
      }
      // Stage 1: mag field propagation for tracks with pstep<safety
      if (stage<=1) {
         // REMAINING ONLY CHARGED TRACKS IN MAG. FIELD
         // Select tracks that undergo safe steps to physics process
         if (fPstepV[itr]<fSafetyV[itr]) {
            fStatusV[itr] = kPhysics;
            step = fPstepV[itr];
            // This category of tracks can make the full step to the physics process
            PropagateInVolumeSingle(itr, step);
            MarkRemoved(itr);
            continue;
         }
      }
      // Stage 2: Remaining tracks trying to cross
      if (stage<=2) {
         c = Curvature(itr);
         if (0.25*c*fSnextV[itr]<1E-6 && fSnextV[itr]<1E-3 && fSnextV[itr]<fPstepV[itr]-1E-6) {
            // Propagate with snext and check if we crossed
            if (fIzeroV[itr]>10) fSnextV[itr] = 1.E-3;
            icrossed += PropagateInFieldSingle(itr, fSnextV[itr]+10*gTolerance, kTRUE);
            if (fSnextV[itr]<1.E-6) fIzeroV[itr]++;
            else fIzeroV[itr] = 0;
            // Crossing tracks have the correct status and are now marked for removal
            gPropagator->fNsnextSteps++;
            if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
            MarkRemoved(itr);
            continue; // -> to next track
         }
         // Track has safety<pstep but next boundary not close enough.
         // We propagate in field with the safety value.
         // Protection in case we have a very small safey after entering a new volume
         const Double_t fs = 0.1;
//         if (fSafetyV[itr]<1.E-6) {
         if (/*fSafetyV[itr] < 0.001 &&*/ fSafetyV[itr] < fs*fSnextV[itr]) {
            // Track getting away from boundary. Work to be done here
            // In principle we need a safety value for the content of the current volume only
            // This does not behave well on corners...
            // ... so we peek a small value and check if this crosses, than recompute safety
//            fSafetyV[itr] = 1.E-3;

            Double_t delta = fs*fSnextV[itr];
            if (c*fSnextV[itr]<0.1) delta = 0.5*fSnextV[itr]*fSnextV[itr]*c;
            else if (c*fSnextV[itr]>10.) delta = TMath::Max(1/c, delta);
            fSafetyV[itr] = TMath::Max(fSafetyV[itr], delta);
            fIzeroV[itr]++;
            if (fIzeroV[itr] > 10) {fSafetyV[itr] = 0.5*fSnextV[itr]; fIzeroV[itr]=0;}
            icrossed += PropagateInFieldSingle(itr, fSafetyV[itr], kTRUE);
         } else {
            if (fIzeroV[itr] > 10) {
               // Propagate with snext
               icrossed += PropagateInFieldSingle(itr, fSnextV[itr]+10*gTolerance, kTRUE);
               fIzeroV[itr] = 0;
               gPropagator->fNsnextSteps++;
               if (fPathV[itr]->IsOutside()) fStatusV[itr] = kExitingSetup;
               continue; // -> to next track
            }
            if (fSafetyV[itr]<1.E-3) fIzeroV[itr]++;
            // Propagate with safety without checking crossing
            icrossed += PropagateInFieldSingle(itr, fSafetyV[itr], kFALSE);
         }
         gPropagator->fNsafeSteps++;
         if (fPathV[itr]->IsOutside()) {
            fStatusV[itr] = kExitingSetup;
            continue;
         }
      }
   }
   // Compact remaining tracks and move the removed oned to the output container
   if (!fCompact) Compact(&output);
//   if (fNtracks) ComputeTransportLength(fNtracks);
   return icrossed;
}

//______________________________________________________________________________
Double_t GeantTrack_v::Curvature(Int_t i) const
{
// Curvature
   if (fChargeV[i]==0) return 0.;
   return TMath::Abs(kB2C*gPropagator->fBmag/Pt(i));
}

//______________________________________________________________________________
Int_t GeantTrack_v::PostponeTracks(GeantTrack_v &output)
{
// Postpone transport of remaining tracks and copy them to the output.
   Int_t npostponed = fNtracks;
   for (Int_t itr=0; itr<fNtracks; itr++) fStatusV[itr] = kPostponed;
   // Move these tracks to the output container
   output.AddTracks(*this, 0, fNtracks-1);
   RemoveTracks(0, fNtracks-1);
   Clear();
   return npostponed;
}

//______________________________________________________________________________
Int_t GeantTrack_v::PostponeTrack(Int_t itr, GeantTrack_v &output)
{
   // Postpone transport of a track and copy it to the output.
   // Returns where in the output the track was added.

   fStatusV[itr] = kPostponed;
   // Move these tracks to the output container
   Int_t new_itr = output.AddTrack(*this, itr);
   MarkRemoved(itr);
   return new_itr;
}
