#include "GeantTrack.h"
#include "globals.h"
#include <execinfo.h>

#if USE_VECGEOM_NAVIGATOR == 1
 #pragma message("Compiling against VecGeom")
 #include "backend/Backend.h"
 #include "navigation/SimpleNavigator.h"
 #include "volumes/PlacedVolume.h" // equivalent of TGeoNode
 #include "base/Vector3D.h"
 #include "base/Transformation3D.h"
 #include "base/Global.h"
 #include "management/GeoManager.h"
 #ifdef CROSSCHECK
  #include "TGeoNavigator.h"
  #include "TGeoNode.h"
 #endif
//#endif
#else // TGeoNavigator as default
 #pragma message("Compiling against TGeo")
 #include <iostream>
 #include "TGeoNavigator.h"
 #include "TGeoNode.h"
#endif
#include "TGeoManager.h"

#include "WorkloadManager.h"

#ifdef __STAT_DEBUG_TRK
#include "GeantTrackStat.h"
#endif

#include "GeantThreadData.h"
//#include "TGeoHelix.h"
#ifdef GEANT_NVCC
#warning "ConstFieldHelixStepper required but not compileable in NVCC."
#else
#include "ConstFieldHelixStepper.h"
#endif
#include "GeantScheduler.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h>
#else
#include "mm_malloc.h"
#endif
#include <cassert>

#ifdef GEANT_CUDA_DEVICE_BUILD
__constant__ double gTolerance;
#else
const Double_t gTolerance = TGeoShape::Tolerance();
#endif

ClassImp(GeantTrack)

    //______________________________________________________________________________
    GeantTrack::GeantTrack()
    : fEvent(-1), fEvslot(-1), fParticle(-1), fPDG(0), fG5code(0), fEindex(0), fCharge(0),
      fProcess(-1), fVindex(0), fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0),
      fYpos(0), fZpos(0), fXdir(0), fYdir(0), fZdir(0), fP(0), fE(0), fTime(0), fEdep(0),
      fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fFrombdr(false), fPending(false), fPath(0),
      fNextpath(0) {
  // Dummy constructor
}

/* Obtain a backtrace and print it to stdout. */
void print_trace(void) {
  void *array[10];
  size_t size;
  char **strings;
  size_t i;

  size = backtrace(array, 10);
  strings = backtrace_symbols(array, size);

  printf("Obtained %zd stack frames.\n", size);

  for (i = 0; i < size; i++)
    printf("%s\n", strings[i]);

  free(strings);
}
//______________________________________________________________________________
GeantTrack::GeantTrack(Int_t ipdg)
    : fEvent(-1), fEvslot(-1), fParticle(-1), fPDG(ipdg), fG5code(0), fEindex(0), fCharge(0),
      fProcess(-1), fVindex(0), fNsteps(0), fSpecies(kHadron), fStatus(kAlive), fMass(0), fXpos(0),
      fYpos(0), fZpos(0), fXdir(0), fYdir(0), fZdir(0), fP(0), fE(0), fTime(0), fEdep(0),
      fPstep(1.E20), fStep(0), fSnext(0), fSafety(0), fFrombdr(false), fPending(false), fPath(0),
      fNextpath(0) {
  // Constructor
  Int_t maxdepth = GeantPropagator::Instance()->fMaxDepth;
  fPath = VolumePath_t::MakeInstance(maxdepth);
  fNextpath = VolumePath_t::MakeInstance(maxdepth);
}

//______________________________________________________________________________
GeantTrack::GeantTrack(const GeantTrack &other)
    : fEvent(other.fEvent), fEvslot(other.fEvslot), fParticle(other.fParticle), fPDG(other.fPDG),
      fG5code(other.fG5code), fEindex(other.fEindex), fCharge(other.fCharge),
      fProcess(other.fProcess), fVindex(other.fVindex), fNsteps(other.fNsteps),
      fSpecies(other.fSpecies), fStatus(other.fStatus), fMass(other.fMass), fXpos(other.fXpos),
      fYpos(other.fYpos), fZpos(other.fZpos), fXdir(other.fXdir), fYdir(other.fYdir),
      fZdir(other.fZdir), fP(other.fP), fE(other.fE), fTime(other.fTime), fEdep(other.fEdep),
      fPstep(other.fPstep), fStep(other.fStep), fSnext(other.fSnext), fSafety(other.fSafety),
      fFrombdr(other.fFrombdr), fPending(other.fPending), fPath(0), fNextpath(0) {
  // Copy constructor
  Int_t maxdepth = GeantPropagator::Instance()->fMaxDepth;
  fPath = VolumePath_t::MakeInstance(maxdepth);
  fNextpath = VolumePath_t::MakeInstance(maxdepth);
  *fPath = *other.fPath;
  *fNextpath = *other.fNextpath;
}

//______________________________________________________________________________
GeantTrack &GeantTrack::operator=(const GeantTrack &other) {
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
    fVindex = other.fVindex;
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
    fTime = other.fTime;
    fEdep = other.fEdep;
    fPstep = other.fPstep;
    fStep = other.fStep;
    fSnext = other.fSnext;
    fSafety = other.fSafety;
    fFrombdr = other.fFrombdr;
    fPending = other.fPending;
    Int_t maxdepth = GeantPropagator::Instance()->fMaxDepth;
    fPath = VolumePath_t::MakeInstance(maxdepth);
    fNextpath = VolumePath_t::MakeInstance(maxdepth);
    *fPath = *other.fPath;
    *fNextpath = *other.fNextpath;
  }
  return *this;
}

//______________________________________________________________________________
GeantTrack::~GeantTrack() {
  // Destructor.
  VolumePath_t::ReleaseInstance(fPath);
  VolumePath_t::ReleaseInstance(fNextpath);
}

//______________________________________________________________________________
void GeantTrack::ReadFromVector(const GeantTrack_v &arr, Int_t i) {
  // Fill track from array
  fEvent = arr.fEventV[i];
  fEvslot = arr.fEvslotV[i];
  fParticle = arr.fParticleV[i];
  fPDG = arr.fPDGV[i];
  fG5code = arr.fG5codeV[i];
  fEindex = arr.fEindexV[i];
  fCharge = arr.fChargeV[i];
  fProcess = arr.fProcessV[i];
  fVindex = arr.fVindexV[i];
  fNsteps = arr.fNstepsV[i];
  fSpecies = arr.fSpeciesV[i];
  fStatus = arr.fStatusV[i];
  fMass = arr.fMassV[i];
  fXpos = arr.fXposV[i];
  fYpos = arr.fYposV[i];
  fZpos = arr.fZposV[i];
  fXdir = arr.fXdirV[i];
  fYdir = arr.fYdirV[i];
  fZdir = arr.fZdirV[i];
  fP = arr.fPV[i];
  fE = arr.fEV[i];
  fTime = arr.fTimeV[i];
  fEdep = arr.fEdepV[i];
  fPstep = arr.fPstepV[i];
  fStep = arr.fStepV[i];
  fSnext = arr.fSnextV[i];
  fSafety = arr.fSafetyV[i];
  fFrombdr = arr.fFrombdrV[i];
  fPending = arr.fPendingV[i];
  //   if (!fPath) fPath = wm->NavStates()->borrow();
  *fPath = *arr.fPathV[i];
  //   if (!fNextpath) fNextpath = wm->NavStates()->borrow();
  *fNextpath = *arr.fNextpathV[i];
}

//______________________________________________________________________________
void GeantTrack::Clear(Option_t *) {
  // Resets track content.
  fEvent = -1;
  fEvslot = -1;
  fParticle = -1;
  fPDG = 0;
  fG5code = 0;
  fEindex = 0;
  fCharge = 0;
  fProcess = -1;
  fVindex = 0;
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
  fTime = 0.;
  fEdep = 0;
  fPstep = 1.E20;
  fStep = 0.;
  fSnext = 0.;
  fSafety = 0.;
  fFrombdr = false;
  fPending = false;
}

//______________________________________________________________________________
Double_t GeantTrack::Curvature() const {
  // Curvature
  if (fCharge == 0)
    return 0.;
  const Double_t tiny = 1.E-30;
  return Math::Abs(kB2C * fCharge * gPropagator->fBmag / (Pt() + tiny));
}

//______________________________________________________________________________
TGeoVolume *GeantTrack::GetVolume() const {
  // Current volume the track is into
  return ((TGeoVolume*)gGeoManager->GetListOfVolumes()->At(fVindex));
}

//______________________________________________________________________________
TGeoVolume *GeantTrack::GetNextVolume() const {
  // Next volume the track is entering
  return fNextpath->GetCurrentNode()->GetVolume();
}

//______________________________________________________________________________
TGeoMaterial *GeantTrack::GetMaterial() const {
  // Current material the track is into
  TGeoMedium *med = GetVolume()->GetMedium();
  if (!med)
    return 0;
  return med->GetMaterial();
}

//______________________________________________________________________________
void GeantTrack::SetPath(VolumePath_t const *const path) {
  // Set path.
  *fPath = *path;
}

//______________________________________________________________________________
void GeantTrack::SetNextPath(VolumePath_t const *const path) {
  // Set next path.
  *fNextpath = *path;
}

//______________________________________________________________________________
void GeantTrack::Print(Int_t) const {
  TString spath;
  //   if (path) path->GetPath(spath);
  Printf("=== Track %d (ev=%d): Process=%d, pstep=%g Charge=%d  Position:(%f,%f,%f) Dir:(%f,%f,%f) "
         "P:%g E:%g snext=%g safety=%g nsteps=%d",
         fParticle, fEvent, fProcess, fPstep, fCharge, fXpos, fYpos, fZpos, fXdir, fYdir, fZdir,
         P(), fE, fSnext, fSafety, fNsteps);
}

ClassImp(GeantTrack_v)

    //______________________________________________________________________________
    GeantTrack_v::GeantTrack_v()
    : fNtracks(0), fMaxtracks(0), fNselected(0), fHoles(0), fSelected(0), fCompact(true),
      fMixed(false), fMaxDepth(0), fBufSize(0), fVPstart(0), fBuf(0), fEventV(0), fEvslotV(0),
      fParticleV(0), fPDGV(0), fG5codeV(0), fEindexV(0), fChargeV(0), fProcessV(0), fVindexV(0),
      fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0),
      fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0), fStepV(0), fSnextV(0),
      fSafetyV(0), fFrombdrV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
// Dummy ctor.
#ifdef __STAT_DEBUG_TRK
  fStat.InitArrays(gPropagator->fNevents);
#endif
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(Int_t size, Int_t maxdepth)
    : fNtracks(0), fMaxtracks(0), fNselected(0), fHoles(0), fSelected(0), fCompact(true),
      fMixed(false), fMaxDepth(maxdepth), fBufSize(0), fVPstart(0), fBuf(0), fEventV(0),
      fEvslotV(0), fParticleV(0), fPDGV(0), fG5codeV(0), fEindexV(0), fChargeV(0), fProcessV(0),
      fVindexV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0),
      fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0),
      fStepV(0), fSnextV(0), fSafetyV(0), fFrombdrV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
// Constructor with maximum capacity.
#ifdef __STAT_DEBUG_TRK
  fStat.InitArrays(gPropagator->fNevents);
#endif
  Resize(size);
}

//______________________________________________________________________________
GeantTrack_v *GeantTrack_v::MakeInstanceAt(void *addr, unsigned int nTracks, Int_t maxdepth)
{
   return new (addr) GeantTrack_v(addr,nTracks,maxdepth);
}


//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(void *addr, unsigned int nTracks, Int_t maxdepth)
    : fNtracks(0), fMaxtracks(nTracks), fNselected(0), fHoles(0), fSelected(0), fCompact(true),
      fMixed(false), fMaxDepth(maxdepth), fBufSize(0), fVPstart(0), fBuf(0), fEventV(0),
      fEvslotV(0), fParticleV(0), fPDGV(0), fG5codeV(0), fEindexV(0), fChargeV(0), fProcessV(0),
      fVindexV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0),
      fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0),
      fStepV(0), fSnextV(0), fSafetyV(0), fFrombdrV(0), fPendingV(0), fPathV(0), fNextpathV(0) {

// Constructor with maximum capacity.
#ifdef __STAT_DEBUG_TRK
  fStat.InitArrays(gPropagator->fNevents);
#endif

  fBuf = ((char*)addr) + sizeof(GeantTrack_v);
  fBufSize = SizeOfInstance(nTracks, maxdepth);
  memset(fBuf, 0, fBufSize);
  AssignInBuffer(fBuf, nTracks);
  memset(fPathV, 0, nTracks * sizeof(VolumePath_t *));
  memset(fNextpathV, 0, nTracks * sizeof(VolumePath_t *));
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(const GeantTrack_v &track_v)
    : fNtracks(0), fMaxtracks(track_v.fMaxtracks), fNselected(track_v.fNselected), fHoles(0),
      fSelected(0), fCompact(track_v.fCompact), fMixed(track_v.fMixed),
      fMaxDepth(track_v.fMaxDepth), fBufSize(track_v.fBufSize), fVPstart(0), fBuf(0), fEventV(0),
      fEvslotV(0), fParticleV(0), fPDGV(0), fG5codeV(0), fEindexV(0), fChargeV(0), fProcessV(0),
      fVindexV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0),
      fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0),
      fStepV(0), fSnextV(0), fSafetyV(0), fFrombdrV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
// Copy constructor
#ifdef __STAT_DEBUG_TRK
  fStat.InitArrays(gPropagator->fNevents);
#endif
#ifndef GEANT_CUDA_DEVICE_BUILD
  fNtracks.store(track_v.fNtracks);
#else
  fNtracks = track_v.fNtracks;
#endif
  fBuf = (char *)_mm_malloc(fBufSize, ALIGN_PADDING);
  memcpy(fBuf, track_v.fBuf, fBufSize);
  AssignInBuffer(&fBuf[0], fMaxtracks);
}

//______________________________________________________________________________
GeantTrack_v &GeantTrack_v::operator=(const GeantTrack_v &track_v) {
  // Assignment operator
  if (&track_v != this) {
#ifndef GEANT_CUDA_DEVICE_BUILD
    fNtracks.store(track_v.fNtracks);
#else
    fNtracks = track_v.fNtracks;
#endif
    Int_t size = track_v.fMaxtracks;
    fMaxDepth = track_v.fMaxDepth;
    fBufSize = track_v.fBufSize;
    if (fMaxtracks < size) {
      _mm_free(fBuf);
      fBuf = (char *)_mm_malloc(fBufSize, ALIGN_PADDING);
    }
    fMaxtracks = size;
    fNselected = track_v.fNselected;
    fHoles = 0;
    fSelected = 0;
    fCompact = track_v.fCompact;
    fMixed = track_v.fMixed;
    memcpy(fBuf, track_v.fBuf, size * sizeof(GeantTrack));
    AssignInBuffer(&fBuf[0], size);
#ifdef __STAT_DEBUG_TRK
    fStat.InitArrays(gPropagator->fNevents);
#endif
  }
  return *this;
}

//______________________________________________________________________________
GeantTrack_v::~GeantTrack_v() {
  // Destructor.
  Int_t ntracks = GetNtracks();
  for (auto i = 0; i < ntracks; ++i) {
    VolumePath_t::ReleaseInstance(fPathV[i]);
    VolumePath_t::ReleaseInstance(fNextpathV[i]);
  }
  _mm_free(fBuf);
}

//______________________________________________________________________________
void GeantTrack_v::AssignInBuffer(char *buff, Int_t size) {
  // Assign all internal class arrays in the supplied buffer, padded by supplied
  // size.
  const Int_t size_intn = size * sizeof(Int_t);
  const Int_t size_doublen = size * sizeof(Double_t);
  const Int_t size_booln = size * sizeof(Bool_t);
  char *buf = buff;
  fEventV = (Int_t *)buf;
  buf += size_intn;
  fEvslotV = (Int_t *)buf;
  buf += size_intn;
  fParticleV = (Int_t *)buf;
  buf += size_intn;
  fPDGV = (Int_t *)buf;
  buf += size_intn;
  fG5codeV = (Int_t *)buf;
  buf += size_intn;
  fEindexV = (Int_t *)buf;
  buf += size_intn;
  fChargeV = (Int_t *)buf;
  buf += size_intn;
  fProcessV = (Int_t *)buf;
  buf += size_intn;
  fVindexV = (Int_t *)buf;
  buf += size_intn;
  fNstepsV = (Int_t *)buf;
  buf += size_intn;
  fSpeciesV = (Species_t *)buf;
  buf += size * sizeof(Species_t);
  fStatusV = (TrackStatus_t *)buf;
  buf += size * sizeof(TrackStatus_t);
  fMassV = (Double_t *)buf;
  buf += size_doublen;
  fXposV = (Double_t *)buf;
  buf += size_doublen;
  fYposV = (Double_t *)buf;
  buf += size_doublen;
  fZposV = (Double_t *)buf;
  buf += size_doublen;
  fXdirV = (Double_t *)buf;
  buf += size_doublen;
  fYdirV = (Double_t *)buf;
  buf += size_doublen;
  fZdirV = (Double_t *)buf;
  buf += size_doublen;
  fPV = (Double_t *)buf;
  buf += size_doublen;
  fEV = (Double_t *)buf;
  buf += size_doublen;
  fTimeV = (Double_t *)buf;
  buf += size_doublen;
  fEdepV = (Double_t *)buf;
  buf += size_doublen;
  fPstepV = (Double_t *)buf;
  buf += size_doublen;
  fStepV = (Double_t *)buf;
  buf += size_doublen;
  fSnextV = (Double_t *)buf;
  buf += size_doublen;
  fSafetyV = (Double_t *)buf;
  buf += size_doublen;
  fFrombdrV = (Bool_t *)buf;
  buf += size_booln;
  fPendingV = (Bool_t *)buf;
  buf += size_booln;
  fPathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  fNextpathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  fVPstart = buf;
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  // Allocate VolumePath_t objects in the reserved buffer space
  for (auto i = 0; i < 2 * size; ++i)
    VolumePath_t::MakeInstanceAt(fMaxDepth, buf + i * size_vpath);
  buf += 2 * size * size_vpath;
  size_t size_bits = BitSet::SizeOfInstance(size);
  fHoles = BitSet::MakeInstanceAt(size, buf);
  buf += size_bits;
  fSelected = BitSet::MakeInstanceAt(size, buf);
}

//______________________________________________________________________________
void GeantTrack_v::CopyToBuffer(char *buff, Int_t size) {
  // Copy existing track arrays into new buffer, padded by supplied size
  Int_t ntracks = GetNtracks();
  const Int_t size_int = ntracks * sizeof(Int_t);
  const Int_t size_double = ntracks * sizeof(Double_t);
  const Int_t size_intn = size * sizeof(Int_t);
  const Int_t size_doublen = size * sizeof(Double_t);
  char *buf = buff;
  memcpy_align(buf, fEventV, size_int);
  fEventV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fEvslotV, size_int);
  fEvslotV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fParticleV, size_int);
  fParticleV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fPDGV, size_int);
  fPDGV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fG5codeV, size_int);
  fG5codeV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fEindexV, size_int);
  fEindexV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fChargeV, size_int);
  fChargeV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fProcessV, size_int);
  fProcessV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fVindexV, size_int);
  fVindexV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fNstepsV, size_int);
  fNstepsV = (Int_t *)buf;
  buf += size_intn;
  memcpy_align(buf, fSpeciesV, ntracks * sizeof(Species_t));
  fSpeciesV = (Species_t *)buf;
  buf += size * sizeof(Species_t);
  memcpy_align(buf, fStatusV, ntracks * sizeof(TrackStatus_t));
  fStatusV = (TrackStatus_t *)buf;
  buf += size * sizeof(TrackStatus_t);
  memcpy_align(buf, fMassV, size_double);
  fMassV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fXposV, size_double);
  fXposV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fYposV, size_double);
  fYposV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fZposV, size_double);
  fZposV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fXdirV, size_double);
  fXdirV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fYdirV, size_double);
  fYdirV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fZdirV, size_double);
  fZdirV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fPV, size_double);
  fPV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fEV, size_double);
  fEV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fTimeV, size_double);
  fTimeV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fEdepV, size_double);
  fEdepV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fPstepV, size_double);
  fPstepV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fStepV, size_double);
  fStepV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fSnextV, size_double);
  fSnextV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fSafetyV, size_double);
  fSafetyV = (Double_t *)buf;
  buf += size_doublen;
  memcpy_align(buf, fFrombdrV, ntracks * sizeof(Bool_t));
  fFrombdrV = (Bool_t *)buf;
  buf += size * sizeof(Bool_t);
  memcpy_align(buf, fPendingV, ntracks * sizeof(Bool_t));
  fPendingV = (Bool_t *)buf;
  buf += size * sizeof(Bool_t);
  //   memcpy_align(buf, fPathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **pathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  //   memcpy_align(buf, fNextpathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **nextpathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  fVPstart = buf;
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  // Allocate VolumePath_t objects in the reserved buffer space
  for (auto i = 0; i < 2 * size; ++i)
    VolumePath_t::MakeInstanceAt(fMaxDepth, fVPstart + i * size_vpath);
  // Copy existing path and nextpath into new buffer
  for (auto i = 0; i < ntracks; ++i) {
    pathV[i] = reinterpret_cast<VolumePath_t *>(fVPstart + i * size_vpath);
    nextpathV[i] = reinterpret_cast<VolumePath_t *>(fVPstart + (size + i) * size_vpath);
    fPathV[i]->CopyTo(pathV[i]);
    fNextpathV[i]->CopyTo(nextpathV[i]);
    VolumePath_t::ReleaseInstance(fPathV[i]);
    VolumePath_t::ReleaseInstance(fNextpathV[i]);
  }
  // Set the new pointers to arrays
  fPathV = pathV;
  fNextpathV = nextpathV;
  buf += 2 * size * size_vpath;
  size_t size_bits = BitSet::SizeOfInstance(size);
  BitSet *holes = BitSet::MakeCopyAt(*fHoles, buf);
  BitSet::ReleaseInstance(fHoles);
  fHoles = holes;
  buf += size_bits;
  BitSet *selected = BitSet::MakeInstanceAt(size, buf);
  BitSet::ReleaseInstance(fSelected);
  fSelected = selected;
}

//______________________________________________________________________________
Bool_t GeantTrack_v::IsSame(const GeantTrack_v &tr1, Int_t i1, const GeantTrack_v &tr2, Int_t i2) {
  // Compare two tracks.
  Long64_t chk1, chk2;
  chk1 = tr1.fEventV[i1] + tr1.fEvslotV[i1] + tr1.fParticleV[i1] + tr1.fPDGV[i1] +
         tr1.fG5codeV[i1] + tr1.fEindexV[i1] + tr1.fChargeV[i1] + tr1.fProcessV[i1] +
         tr1.fVindexV[i1] + tr1.fNstepsV[i1] + (Long64_t)tr1.fSpeciesV[i1] +
         (Long64_t)tr1.fStatusV[i1];
  chk2 = tr2.fEventV[i2] + tr2.fEvslotV[i2] + tr2.fParticleV[i2] + tr2.fPDGV[i2] +
         tr2.fG5codeV[i2] + tr2.fEindexV[i2] + tr2.fChargeV[i2] + tr2.fProcessV[i2] +
         tr2.fVindexV[i2] + tr2.fNstepsV[i2] + (Long64_t)tr2.fSpeciesV[i2] +
         (Long64_t)tr2.fStatusV[i2];
  if (chk1 != chk2)
    return false;
  Double_t dchk1, dchk2;
  dchk1 = (Long64_t)tr1.fMassV[i1] + tr1.fXposV[i1] + tr1.fYposV[i1] + tr1.fZposV[i1] +
          tr1.fXdirV[i1] + tr1.fYdirV[i1] + tr1.fZdirV[i1] + tr1.fPV[i1] + tr1.fEdepV[i1] +
          tr1.fEV[i1] + tr1.fPstepV[i1] + tr1.fStepV[i1] + tr1.fSnextV[i1] + tr1.fSafetyV[i1];
  dchk2 = (Long64_t)tr2.fMassV[i2] + tr2.fXposV[i2] + tr2.fYposV[i2] + tr2.fZposV[i2] +
          tr2.fXdirV[i2] + tr2.fYdirV[i2] + tr2.fZdirV[i2] + tr2.fPV[i2] + tr2.fEdepV[i2] +
          tr2.fEV[i2] + tr2.fPstepV[i2] + tr2.fStepV[i2] + tr2.fSnextV[i2] + tr2.fSafetyV[i2];
  if (!TMath::AreEqualAbs(dchk1, dchk2, 1.E-10))
    return false;
  if (tr1.fPendingV[i1] != tr2.fPendingV[i2])
    return false;
  return true;
}

//______________________________________________________________________________
void GeantTrack_v::CheckTracks() {
  //  for (Int_t i=0; i<fMaxtracks; ++i)
  //     if (fNextpathV[i]) fNextpathV[i]->SetClient(gPropagator);
}

//______________________________________________________________________________
size_t GeantTrack_v::SizeOfInstance(size_t nTracks, size_t maxdepth) {
   // return the contiguous memory size needed to hold a GeantTrack_v

   size_t size = round_up_align(nTracks);
   size_t size_nav = 2 * size * VolumePath_t::SizeOfInstance(maxdepth);
   size_t size_bits = 2 * BitSet::SizeOfInstance(size);

   return size * sizeof(GeantTrack) + size_nav + size_bits;
}

//______________________________________________________________________________
void GeantTrack_v::Resize(Int_t newsize) {
  // Resize the container.
  Int_t size = round_up_align(newsize);
  if (size < GetNtracks()) {
    Printf("Error: Cannot resize to less than current track content");
    return;
  }
  fBufSize = SizeOfInstance(size, fMaxDepth);
  if (!fCompact)
    Compact();

  char *buf = (char *)_mm_malloc(fBufSize, ALIGN_PADDING);
  memset(buf, 0, fBufSize);
  fMaxtracks = size;
  if (!fBuf) {
    // All arrays are contiguous in a single buffer and aligned with the
    // same padding ALIGN_PADDING
    fBuf = buf;
    AssignInBuffer(buf, size);
    memset(fPathV, 0, size * sizeof(VolumePath_t *));
    memset(fNextpathV, 0, size * sizeof(VolumePath_t *));
  } else {
    // Resize container
    CopyToBuffer(buf, size);
    _mm_free(fBuf);
    fBuf = buf;
  }
  fHoles->ResetAllBits();
  fSelected->ResetAllBits();
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Int_t GeantTrack_v::AddTrack(GeantTrack &track, Bool_t /*import*/) {
  // Add new track to the array. If addition is done on top of non-compact array,
  // the track will be inserted without updating the number of tracks. If track is
  // imported just copy the pointers to the navigation states and reset the sources.
  // Returns the location where the track was added.
  Int_t itrack = GetNtracks();
  if (!fCompact)
    itrack = fHoles->FirstSetBit();
  if (itrack == fMaxtracks) {
#ifndef GEANT_CUDA_DEVICE_BUILD
    Resize(2 * fMaxtracks);
#else
    printf("Error in GeantTrack_v::AddTrack, resizing is not supported in device code\n");
#endif
  }
  fHoles->ResetBitNumber(itrack);
  fSelected->ResetBitNumber(itrack);
  fEventV[itrack] = track.fEvent;
  fEvslotV[itrack] = track.fEvslot;
  fParticleV[itrack] = track.fParticle;
  fPDGV[itrack] = track.fPDG;
  fG5codeV[itrack] = track.fG5code;
  fEindexV[itrack] = track.fEindex;
  fChargeV[itrack] = track.fCharge;
  fProcessV[itrack] = track.fProcess;
  fVindexV[itrack] = track.fVindex;
  fNstepsV[itrack] = track.fNsteps;
  fSpeciesV[itrack] = track.fSpecies;
  fStatusV[itrack] = track.fStatus;
  fMassV[itrack] = track.fMass;
  fXposV[itrack] = track.fXpos;
  fYposV[itrack] = track.fYpos;
  fZposV[itrack] = track.fZpos;
  fXdirV[itrack] = track.fXdir;
  fYdirV[itrack] = track.fYdir;
  fZdirV[itrack] = track.fZdir;
  fPV[itrack] = track.fP;
  fEV[itrack] = track.fE;
  fTimeV[itrack] = track.fTime;
  fEdepV[itrack] = track.fEdep;
  fPstepV[itrack] = track.fPstep;
  fStepV[itrack] = track.fStep;
  fSnextV[itrack] = track.fSnext;
  fSafetyV[itrack] = track.fSafety;
  fFrombdrV[itrack] = track.fFrombdr;
  fPendingV[itrack] = track.fPending;
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  track.fPath->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] =
      reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  track.fNextpath->CopyTo(fNextpathV[itrack]);
  fNtracks++;
#ifdef __STAT_DEBUG_TRK
  fStat.fNtracks[fEvslotV[itrack]]++;
#endif
  return itrack;
}

//______________________________________________________________________________
Int_t GeantTrack_v::AddTrackSync(GeantTrack &track) {
  // Add track in a concurrent way. Assumes that this array
  // Is currently being filled while held by the basket manager and NOT being
  // transported.
  // The array has to be compact and should have enough alocated space.
  // Returns the location where the track was added.
  assert(fCompact);
  assert(GetNtracks() < fMaxtracks);
  Int_t itrack = fNtracks++;
  fEventV[itrack] = track.fEvent;
  fEvslotV[itrack] = track.fEvslot;
  fParticleV[itrack] = track.fParticle;
  fPDGV[itrack] = track.fPDG;
  fG5codeV[itrack] = track.fG5code;
  fEindexV[itrack] = track.fEindex;
  fChargeV[itrack] = track.fCharge;
  fProcessV[itrack] = track.fProcess;
  fVindexV[itrack] = track.fVindex;
  fNstepsV[itrack] = track.fNsteps;
  fSpeciesV[itrack] = track.fSpecies;
  fStatusV[itrack] = track.fStatus;
  fMassV[itrack] = track.fMass;
  fXposV[itrack] = track.fXpos;
  fYposV[itrack] = track.fYpos;
  fZposV[itrack] = track.fZpos;
  fXdirV[itrack] = track.fXdir;
  fYdirV[itrack] = track.fYdir;
  fZdirV[itrack] = track.fZdir;
  fPV[itrack] = track.fP;
  fEV[itrack] = track.fE;
  fTimeV[itrack] = track.fTime;
  fEdepV[itrack] = track.fEdep;
  fPstepV[itrack] = track.fPstep;
  fStepV[itrack] = track.fStep;
  fSnextV[itrack] = track.fSnext;
  fSafetyV[itrack] = track.fSafety;
  fFrombdrV[itrack] = track.fFrombdr;
  fPendingV[itrack] = track.fPending;
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  track.fPath->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] =
      reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  track.fNextpath->CopyTo(fNextpathV[itrack]);
#ifdef __STAT_DEBUG_TRK
  fStat.fNtracks[fEvslotV[itrack]]++;
#endif
  return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::GetTrack(Int_t i, GeantTrack &track) const {
  // Extract a single track from array.
  track.ReadFromVector(*this, i);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Int_t GeantTrack_v::AddTrack(GeantTrack_v &arr, Int_t i, Bool_t /*import*/) {
// Add track from different array
// If addition is done on top of non-compact array,
// the track will be inserted without updating the number of tracks.
// Returns the location where the track was added.
#ifdef VERBOSE
  arr.PrintTrack(i);
#endif
  Int_t itrack = GetNtracks();
  if (!fCompact)
    itrack = fHoles->FirstSetBit();
  if (itrack == fMaxtracks) {
#ifndef GEANT_CUDA_DEVICE_BUILD
    Resize(2 * fMaxtracks);
#else
    printf("Error in GeantTrack_v::AddTrack, resizing is not supported in device code\n");
#endif
  }
  fHoles->ResetBitNumber(itrack);
  fSelected->ResetBitNumber(itrack);

  fEventV[itrack] = arr.fEventV[i];
  fEvslotV[itrack] = arr.fEvslotV[i];
  fParticleV[itrack] = arr.fParticleV[i];
  fPDGV[itrack] = arr.fPDGV[i];
  fG5codeV[itrack] = arr.fG5codeV[i];
  fEindexV[itrack] = arr.fEindexV[i];
  fChargeV[itrack] = arr.fChargeV[i];
  fProcessV[itrack] = arr.fProcessV[i];
  fVindexV[itrack] = arr.fVindexV[i];
  fNstepsV[itrack] = arr.fNstepsV[i];
  fSpeciesV[itrack] = arr.fSpeciesV[i];
  fStatusV[itrack] = arr.fStatusV[i];
  fMassV[itrack] = arr.fMassV[i];
  fXposV[itrack] = arr.fXposV[i];
  fYposV[itrack] = arr.fYposV[i];
  fZposV[itrack] = arr.fZposV[i];
  fXdirV[itrack] = arr.fXdirV[i];
  fYdirV[itrack] = arr.fYdirV[i];
  fZdirV[itrack] = arr.fZdirV[i];
  fPV[itrack] = arr.fPV[i];
  fEV[itrack] = arr.fEV[i];
  fTimeV[itrack] = arr.fTimeV[i];
  fEdepV[itrack] = arr.fEdepV[i];
  fPstepV[itrack] = arr.fPstepV[i];
  fStepV[itrack] = arr.fStepV[i];
  fSnextV[itrack] = arr.fSnextV[i];
  fSafetyV[itrack] = arr.fSafetyV[i];
  fFrombdrV[itrack] = arr.fFrombdrV[i];
  fPendingV[itrack] = arr.fPendingV[i];
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  arr.fPathV[i]->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] =
      reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  arr.fNextpathV[i]->CopyTo(fNextpathV[itrack]);
  fNtracks++;
#ifdef __STAT_DEBUG_TRK
  fStat.fNtracks[arr.fEvslotV[i]]++;
#endif

  return itrack;
}

//______________________________________________________________________________
Int_t GeantTrack_v::AddTrackSync(GeantTrack_v &arr, Int_t i) {
  // Add track from different array in a concurrent way. Assumes that this array
  // Is currently being filled while held by the basket manager and NOT being
  // transported.
  // The array has to be compact and should have enough alocated space.
  // Returns the location where the track was added.
  assert(fCompact);
  assert(GetNtracks() < fMaxtracks);
#ifdef VERBOSE
  arr.PrintTrack(i);
#endif
  // WorkloadManager *wm = WorkloadManager::Instance();
  Int_t itrack = fNtracks++;

  fEventV[itrack] = arr.fEventV[i];
  fEvslotV[itrack] = arr.fEvslotV[i];
  fParticleV[itrack] = arr.fParticleV[i];
  fPDGV[itrack] = arr.fPDGV[i];
  fG5codeV[itrack] = arr.fG5codeV[i];
  fEindexV[itrack] = arr.fEindexV[i];
  fChargeV[itrack] = arr.fChargeV[i];
  fProcessV[itrack] = arr.fProcessV[i];
  fVindexV[itrack] = arr.fVindexV[i];
  fNstepsV[itrack] = arr.fNstepsV[i];
  fSpeciesV[itrack] = arr.fSpeciesV[i];
  fStatusV[itrack] = arr.fStatusV[i];
  fMassV[itrack] = arr.fMassV[i];
  fXposV[itrack] = arr.fXposV[i];
  fYposV[itrack] = arr.fYposV[i];
  fZposV[itrack] = arr.fZposV[i];
  fXdirV[itrack] = arr.fXdirV[i];
  fYdirV[itrack] = arr.fYdirV[i];
  fZdirV[itrack] = arr.fZdirV[i];
  fPV[itrack] = arr.fPV[i];
  fEV[itrack] = arr.fEV[i];
  fTimeV[itrack] = arr.fTimeV[i];
  fEdepV[itrack] = arr.fEdepV[i];
  fPstepV[itrack] = arr.fPstepV[i];
  fStepV[itrack] = arr.fStepV[i];
  fSnextV[itrack] = arr.fSnextV[i];
  fSafetyV[itrack] = arr.fSafetyV[i];
  fFrombdrV[itrack] = arr.fFrombdrV[i];
  fPendingV[itrack] = arr.fPendingV[i];
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  arr.fPathV[i]->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] =
      reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  arr.fNextpathV[i]->CopyTo(fNextpathV[itrack]);
#ifdef __STAT_DEBUG_TRK
  fStat.fNtracks[arr.fEvslotV[i]]++;
#endif
  return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::AddTracks(GeantTrack_v &arr, Int_t istart, Int_t iend, Bool_t /*import*/) {
// Add tracks from different array. Single thread at a time.
#ifdef __STAT_DEBUG_TRK
  for (Int_t i = istart; i <= iend; i++)
    fStat.fNtracks[arr.fEvslotV[i]]++;
#endif
  Int_t ncpy = iend - istart + 1;
  Int_t ntracks = GetNtracks();
  if (ntracks + ncpy >= fMaxtracks) {
    Resize(Math::Max(2 * fMaxtracks, ntracks + ncpy));
  }
  memcpy_align(&fEventV[ntracks], &arr.fEventV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fEvslotV[ntracks], &arr.fEvslotV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fParticleV[ntracks], &arr.fParticleV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fPDGV[ntracks], &arr.fPDGV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fG5codeV[ntracks], &arr.fG5codeV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fEindexV[ntracks], &arr.fEindexV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fChargeV[ntracks], &arr.fChargeV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fProcessV[ntracks], &arr.fProcessV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fVindexV[ntracks], &arr.fVindexV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fNstepsV[ntracks], &arr.fNstepsV[istart], ncpy * sizeof(Int_t));
  memcpy_align(&fSpeciesV[ntracks], &arr.fSpeciesV[istart], ncpy * sizeof(Species_t));
  memcpy_align(&fStatusV[ntracks], &arr.fStatusV[istart], ncpy * sizeof(TrackStatus_t));
  memcpy_align(&fMassV[ntracks], &arr.fMassV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fXposV[ntracks], &arr.fXposV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fYposV[ntracks], &arr.fYposV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fZposV[ntracks], &arr.fZposV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fXdirV[ntracks], &arr.fXdirV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fYdirV[ntracks], &arr.fYdirV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fZdirV[ntracks], &arr.fZdirV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fPV[ntracks], &arr.fPV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fEV[ntracks], &arr.fEV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fTimeV[ntracks], &arr.fTimeV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fEdepV[ntracks], &arr.fEdepV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fPstepV[ntracks], &arr.fPstepV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fStepV[ntracks], &arr.fStepV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fSnextV[ntracks], &arr.fSnextV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fSafetyV[ntracks], &arr.fSafetyV[istart], ncpy * sizeof(Double_t));
  memcpy_align(&fFrombdrV[ntracks], &arr.fFrombdrV[istart], ncpy * sizeof(Bool_t));
  memcpy_align(&fPendingV[ntracks], &arr.fPendingV[istart], ncpy * sizeof(Bool_t));

  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  for (auto i = ntracks, j = istart; i < (ntracks + ncpy); ++i, ++j) {
    fPathV[i] = reinterpret_cast<VolumePath_t *>(fVPstart + i * size_vpath);
    fNextpathV[i] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + i) * size_vpath);
    arr.fPathV[j]->CopyTo(fPathV[i]);
    arr.fNextpathV[j]->CopyTo(fNextpathV[i]);
  }
  fSelected->ResetBitNumber(ntracks + ncpy - 1);
  fHoles->ResetBitNumber(ntracks + ncpy - 1);
  fNtracks += ncpy;
}

//______________________________________________________________________________
void GeantTrack_v::SwapTracks(Int_t i, Int_t j) {
  // Swap two tracks in the container
  Double_t tdbl;
  Int_t tint;
  Bool_t tbool;
  VolumePath_t *tptr;
  tint = fEventV[i];
  fEventV[i] = fEventV[j];
  fEventV[j] = tint;
  tint = fEvslotV[i];
  fEvslotV[i] = fEvslotV[j];
  fEvslotV[j] = tint;
  tint = fParticleV[i];
  fParticleV[i] = fParticleV[j];
  fParticleV[j] = tint;
  tint = fPDGV[i];
  fPDGV[i] = fPDGV[j];
  fPDGV[j] = tint;
  tint = fG5codeV[i];
  fG5codeV[i] = fG5codeV[j];
  fG5codeV[j] = tint;
  tint = fEindexV[i];
  fEindexV[i] = fEindexV[j];
  fEindexV[j] = tint;
  tint = fChargeV[i];
  fChargeV[i] = fChargeV[j];
  fChargeV[j] = tint;
  tint = fProcessV[i];
  fProcessV[i] = fProcessV[j];
  fProcessV[j] = tint;
  tint = fVindexV[i];
  fVindexV[i] = fVindexV[j];
  fVindexV[j] = tint;
  tint = fNstepsV[i];
  fNstepsV[i] = fNstepsV[j];
  fNstepsV[j] = tint;
  Species_t tspec = fSpeciesV[i];
  fSpeciesV[i] = fSpeciesV[j];
  fSpeciesV[j] = tspec;
  TrackStatus_t tstat = fStatusV[i];
  fStatusV[i] = fStatusV[j];
  fStatusV[j] = tstat;
  tdbl = fMassV[i];
  fMassV[i] = fMassV[j];
  fMassV[j] = tdbl;
  tdbl = fXposV[i];
  fXposV[i] = fXposV[j];
  fXposV[j] = tdbl;
  tdbl = fYposV[i];
  fYposV[i] = fYposV[j];
  fYposV[j] = tdbl;
  tdbl = fZposV[i];
  fZposV[i] = fZposV[j];
  fZposV[j] = tdbl;
  tdbl = fXdirV[i];
  fXdirV[i] = fXdirV[j];
  fXdirV[j] = tdbl;
  tdbl = fYdirV[i];
  fYdirV[i] = fYdirV[j];
  fYdirV[j] = tdbl;
  tdbl = fZdirV[i];
  fZdirV[i] = fZdirV[j];
  fZdirV[j] = tdbl;
  tdbl = fPV[i];
  fPV[i] = fPV[j];
  fPV[j] = tdbl;
  tdbl = fEV[i];
  fEV[i] = fEV[j];
  fEV[j] = tdbl;
  tdbl = fTimeV[i];
  fTimeV[i] = fTimeV[j];
  fTimeV[j] = tdbl;
  tdbl = fEdepV[i];
  fEdepV[i] = fEdepV[j];
  fEdepV[j] = tdbl;
  tdbl = fPstepV[i];
  fPstepV[i] = fPstepV[j];
  fPstepV[j] = tdbl;
  tdbl = fStepV[i];
  fStepV[i] = fStepV[j];
  fStepV[j] = tdbl;
  tdbl = fSnextV[i];
  fSnextV[i] = fSnextV[j];
  fSnextV[j] = tdbl;
  tdbl = fSafetyV[i];
  fSafetyV[i] = fSafetyV[j];
  fSafetyV[j] = tdbl;
  tbool = fFrombdrV[i];
  fFrombdrV[i] = fFrombdrV[j];
  fFrombdrV[j] = tbool;
  tbool = fPendingV[i];
  fPendingV[i] = fPendingV[j];
  fPendingV[j] = tbool;
  tptr = fPathV[i];
  fPathV[i] = fPathV[j];
  fPathV[j] = tptr;
  tptr = fNextpathV[i];
  fNextpathV[i] = fNextpathV[j];
  fNextpathV[j] = tptr;
  Bool_t sel = fSelected->TestBitNumber(j);
  fSelected->SetBitNumber(j, fSelected->TestBitNumber(i));
  fSelected->SetBitNumber(i, sel);
}

//______________________________________________________________________________
void GeantTrack_v::ReplaceTrack(Int_t i, Int_t j) {
  // Replace content of track i with the one of track j
  // WorkloadManager *wm = WorkloadManager::Instance();
  fEventV[i] = fEventV[j];
  fEvslotV[i] = fEvslotV[j];
  fParticleV[i] = fParticleV[j];
  fPDGV[i] = fPDGV[j];
  fG5codeV[i] = fG5codeV[j];
  fEindexV[i] = fEindexV[j];
  fChargeV[i] = fChargeV[j];
  fProcessV[i] = fProcessV[j];
  fVindexV[i] = fVindexV[j];
  fNstepsV[i] = fNstepsV[j];
  fSpeciesV[i] = fSpeciesV[j];
  fStatusV[i] = fStatusV[j];
  fMassV[i] = fMassV[j];
  fXposV[i] = fXposV[j];
  fYposV[i] = fYposV[j];
  fZposV[i] = fZposV[j];
  fXdirV[i] = fXdirV[j];
  fYdirV[i] = fYdirV[j];
  fZdirV[i] = fZdirV[j];
  fPV[i] = fPV[j];
  fEV[i] = fEV[j];
  fTimeV[i] = fTimeV[j];
  fEdepV[i] = fEdepV[j];
  fPstepV[i] = fPstepV[j];
  fStepV[i] = fStepV[j];
  fSnextV[i] = fSnextV[j];
  fSafetyV[i] = fSafetyV[j];
  fFrombdrV[i] = fFrombdrV[j];
  fPendingV[i] = fPendingV[j];
  //   if (!fPathV[i]) fPathV[i] = wm->NavStates()->Borrow();
  //   if (!fNextpathV[i]) fNextpathV[i] = wm->NavStates()->Borrow();
  fPathV[i] = fPathV[j]; // fPathV[j] = 0;
  fNextpathV[i] = fNextpathV[j]; // fNextpathV[j] = 0;
  fSelected->SetBitNumber(i, fSelected->TestBitNumber(j));
}

//______________________________________________________________________________
void GeantTrack_v::DeleteTrack(Int_t /*itr*/) {
  // Delete branch arrays for this track. The track should not have a copy, this has
  // to be called after a killed track is removed by the scheduler.
  //   WorkloadManager *wm = WorkloadManager::Instance();
  //   wm->NavStates()->release(fPathV[itr]);
  //   fPathV[itr] = 0;
  //   wm->NavStates()->release(fNextpathV[itr]);
  //   fNextpathV[itr] = 0;
}

//______________________________________________________________________________
void GeantTrack_v::RemoveTracks(Int_t from, Int_t to) {
// Remove tracks from the container. The method assumes that the tracks were
// copied to another container beforehand.
#ifdef __STAT_DEBUG_TRK
  for (Int_t i = from; i <= to; i++)
    fStat.fNtracks[fEvslotV[i]]--;
#endif
#ifndef GEANT_CUDA_DEVICE_BUILD
  if (!fCompact)
    Printf("RemoveTracks: Not compact");
#endif
  Int_t ntracks = GetNtracks();
  if (to >= ntracks - 1) {
    Int_t nzero = ntracks - from;
    memset(&fPathV[from], 0, nzero * sizeof(VolumePath_t *));
    memset(&fNextpathV[from], 0, nzero * sizeof(VolumePath_t *));
  }
  Int_t ncpy = fNtracks - to - 1;
  memmove(&fEventV[from], &fEventV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fEvslotV[from], &fEvslotV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fParticleV[from], &fParticleV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fPDGV[from], &fPDGV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fG5codeV[from], &fG5codeV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fEindexV[from], &fEindexV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fChargeV[from], &fChargeV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fProcessV[from], &fProcessV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fVindexV[from], &fVindexV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fNstepsV[from], &fNstepsV[to + 1], ncpy * sizeof(Int_t));
  memmove(&fSpeciesV[from], &fSpeciesV[to + 1], ncpy * sizeof(Species_t));
  memmove(&fStatusV[from], &fStatusV[to + 1], ncpy * sizeof(TrackStatus_t));
  memmove(&fMassV[from], &fMassV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fXposV[from], &fXposV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fYposV[from], &fYposV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fZposV[from], &fZposV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fXdirV[from], &fXdirV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fYdirV[from], &fYdirV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fZdirV[from], &fZdirV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fPV[from], &fPV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fEV[from], &fEV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fTimeV[from], &fTimeV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fEdepV[from], &fEdepV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fPstepV[from], &fPstepV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fStepV[from], &fStepV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fSnextV[from], &fSnextV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fSafetyV[from], &fSafetyV[to + 1], ncpy * sizeof(Double_t));
  memmove(&fFrombdrV[from], &fFrombdrV[to + 1], ncpy * sizeof(Bool_t));
  memmove(&fPendingV[from], &fPendingV[to + 1], ncpy * sizeof(Bool_t));
  memmove(&fPathV[from], &fPathV[to + 1], ncpy * sizeof(VolumePath_t *));
  memmove(&fNextpathV[from], &fNextpathV[to + 1], ncpy * sizeof(VolumePath_t *));
  fNtracks -= to - from + 1;
  fSelected->ResetAllBits();
  fNselected = 0;
}

//______________________________________________________________________________
Int_t GeantTrack_v::Compact(GeantTrack_v *moveto) {
  // Compact the holes in the array. Return number of active elements. This will
  // lose the track fields in the holes, so information from the holes has to be
  // copied beforehand
  Int_t ntracks = GetNtracks();
  if (ntracks == 0 || fCompact)
    return 0;
  fCompact = kTRUE;
  Int_t firsthole = fHoles->FirstSetBit();
  while (firsthole < ntracks) {
    Int_t lastactive = fHoles->LastNullBit(ntracks - 1);
    if (lastactive < ntracks) {
      // move last holes (if any)
      if (moveto && (ntracks - lastactive - 1 > 0))
        moveto->AddTracks(*this, lastactive + 1, ntracks - 1, kTRUE);
      ntracks = lastactive + 1;
      if (firsthole == ntracks) {
        SetNtracks(ntracks);
        return ntracks;
      }
    } else {
      // No active tracks left. First copy the hole track to the output
      if (moveto)
        moveto->AddTracks(*this, firsthole, firsthole + ntracks - 1, kTRUE);
      SetNtracks(0);
      return 0;
    }
    // replace content of first hole with the last active track
    if (moveto)
      moveto->AddTrack(*this, firsthole, kTRUE);
    ReplaceTrack(firsthole, lastactive);
    fHoles->SetBitNumber(firsthole, false);
    fHoles->SetBitNumber(lastactive, true);
    firsthole = fHoles->FirstSetBit(firsthole + 1);
    ntracks--;
  }
  fSelected->ResetAllBits();
  fNselected = 0;
  SetNtracks(ntracks);
  return ntracks;
}

//______________________________________________________________________________
Int_t GeantTrack_v::Reshuffle() {
  // Reshuffle tracks according the selection mask. The selected tracks will be
  // moved at the beginning of the array. Tracks should be compacted before.
  if (GetNtracks() == 0)
    return 0;
  fNselected = GetNtracks();
  Int_t firsthole = fSelected->FirstNullBit();
  while (firsthole < fNselected) {
    Int_t lastsel = fSelected->LastSetBit(fNselected - 1);
    if (lastsel >= fNselected)
      return 0;
    fNselected = lastsel + 1;
    if (firsthole == fNselected)
      return fNselected;
    // exchange tracks pointed by firsthole and lastactive
    SwapTracks(firsthole, lastsel);
    fSelected->SetBitNumber(firsthole, true);
    fSelected->SetBitNumber(lastsel, false);
    firsthole = fSelected->FirstNullBit(firsthole + 1);
    fNselected--;
  }
  return fNselected;
}

//______________________________________________________________________________
Bool_t GeantTrack_v::Contains(Int_t evstart, Int_t nevents) const {
  // Check if the array contains tracks from a given event range
  Int_t evend = evstart + nevents;
  Int_t ntracks = GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (fEventV[itr] >= evstart && fEventV[itr] < evend)
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
void GeantTrack_v::Clear(Option_t *) {
  // Clear track content and selections
  fNselected = 0;
  Int_t ntracks = GetNtracks();
  if (ntracks) {
    memset(fPathV, 0, ntracks * sizeof(VolumePath_t *));
    memset(fNextpathV, 0, ntracks * sizeof(VolumePath_t *));
  }
  fHoles->ResetAllBits();
  fSelected->ResetAllBits();
  fCompact = kTRUE;
  SetNtracks(0);
#ifdef __STAT_DEBUG_TRK
  fStat.Reset();
#endif
}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateStraight(Int_t ntracks, Double_t *crtstep) {
  // Propagate first ntracks along a straight line (neutral particles, no mag.
  // field or for last tiny step). The array has to be reshuffled after selecting
  // the interesting particles using Select method.
  // The crossing tracks get masked as holes in the array.

  // Find next volume
  Int_t icrossed = 0;
  for (Int_t i = 0; i < ntracks; i++) {
    if (fFrombdrV[i]) {
      //*fPathV[i] = *fNextpathV[i];
      fStatusV[i] = kBoundary;
      icrossed++;
    }
  }
  for (Int_t i = 0; i < ntracks; i++) {
    fPstepV[i] -= crtstep[i];
    fSafetyV[i] = 0;
    // Change path to reflect the physical volume for the current track; The
    // relevant path is fPath[i] if the frombdr flag is not set or fNextpath[i]
    // otherwise
    fXposV[i] += crtstep[i] * fXdirV[i];
    fYposV[i] += crtstep[i] * fYdirV[i];
    fZposV[i] += crtstep[i] * fZdirV[i];
    fNstepsV[i]++;
#ifdef USE_VECGEOM_NAVIGATOR
//      CheckLocationPathConsistency(i);
#endif
  }
  return icrossed;
}

//______________________________________________________________________________
void GeantTrack_v::PropagateInVolume(Int_t ntracks, const Double_t *crtstep, Int_t tid) {
  // Propagate the selected tracks with crtstep values. The method is to be called
  // only with  charged tracks in magnetic field. The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)
  for (Int_t i = 0; i < ntracks; i++) {
    PropagateInVolumeSingle(i, crtstep[i], tid);
  }
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrack_v::PropagateInVolumeSingle(Int_t i, Double_t crtstep, Int_t /*tid*/) {
  // Propagate the selected track with crtstep value. The method is to be called
  // only with  charged tracks in magnetic field.The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)
  //   Double_t c = 0.;
  //   const Double_t *point = 0;
  //   const Double_t *newdir = 0;
  //   GeantThreadData *td = gPropagator->fThreadData[tid];
  //   TGeoHelix *fieldp = td->fFieldPropagator;
  // Reset relevant variables
  fStatusV[i] = kInFlight;
  fPstepV[i] -= crtstep;
  if (fPstepV[i] < 1.E-10) {
    fPstepV[i] = 0;
    fStatusV[i] = kPhysics;
#ifndef GEANT_CUDA_DEVICE_BUILD
 //   gPropagator->fNphysSteps++;
#endif
  }
  fSafetyV[i] -= crtstep;
  if (fSafetyV[i] < 1.E-10)
    fSafetyV[i] = 0;
#ifndef GEANT_CUDA_DEVICE_BUILD
  else if (!fFrombdrV[i])
      gPropagator->fNsafeSteps++;
#endif
  fSnextV[i] -= crtstep;
  if (fSnextV[i] < 1.E-10) {
    fSnextV[i] = 0;
    if (fFrombdrV[i]) {
      fStatusV[i] = kBoundary;
#ifndef GEANT_CUDA_DEVICE_BUILD
      gPropagator->fNsnextSteps++;
#endif
    }
  }
  fStepV[i] += crtstep;
// Set curvature, charge
//   c = Curvature(i);
//    Double_t dir[3] = {0};
// NOTE: vectorized treatment in TGeoHelix -> todo
//    fieldp->SetXYcurvature(c);
//    fieldp->SetCharge(fChargeV[i]);
//    fieldp->SetHelixStep(std::fabs(Math::TwoPi()*Pz(i)/(c*Pt(i))));
//    fieldp->InitPoint(fXposV[i],fYposV[i],fZposV[i]);
//    fieldp->InitDirection(fXdirV[i],fYdirV[i], fZdirV[i]);
//    fieldp->UpdateHelix();
//    fieldp->Step(crtstep);
//    point = fieldp->GetCurrentPoint();
//    newdir = fieldp->GetCurrentDirection();
//    memcpy(dir,newdir,3*sizeof(Double_t));
//    Math::Normalize(dir);
//    fXposV[i] = point[0]; fYposV[i] = point[1]; fZposV[i] = point[2];
//    fXdirV[i] = dir[0]; fYdirV[i] = dir[1]; fZdirV[i] = dir[2];
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
  // alternative code with lean stepper would be:
  // ( stepper header has to be included )
#ifndef GEANT_NVCC
  geantv::ConstBzFieldHelixStepper stepper(gPropagator->fBmag);
  double posnew[3];
  double dirnew[3];
  stepper.DoStep(fXposV[i], fYposV[i], fZposV[i], fXdirV[i], fYdirV[i], fZdirV[i], fChargeV[i],
                 fPV[i], crtstep, posnew[0], posnew[1], posnew[2], dirnew[0], dirnew[1], dirnew[2]);

  //  maybe normalize direction here
  Math::Normalize(dirnew);
  //      double diffpos =
  //      (xnew-point[0])*(xnew-point[0])+(ynew-point[1])*(ynew-point[1])+(znew-point[2])*(znew-point[2]);
  //      if (diffpos>1.E-4) {
  //         Printf("difference in pos = %g", diffpos);
  //      }
  fXposV[i] = posnew[0];
  fYposV[i] = posnew[1];
  fZposV[i] = posnew[2];
  fXdirV[i] = dirnew[0];
  fYdirV[i] = dirnew[1];
  fZdirV[i] = dirnew[2];
#endif
}

#ifdef USE_VECGEOM_NAVIGATOR
void GeantTrack_v::CheckLocationPathConsistency(Int_t itr) const {
  VECGEOM_NAMESPACE::NavigationState *a = VECGEOM_NAMESPACE::NavigationState::MakeInstance(
      VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
  a->Clear();
  VECGEOM_NAMESPACE::SimpleNavigator nav;
  nav.LocatePoint(VECGEOM_NAMESPACE::GeoManager::Instance().GetWorld(),
                  VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision>(
                      fXposV[itr], fYposV[itr], fZposV[itr]),
                  *a, true);
  if (a->Top() != NULL && a->Top() != fPathV[itr]->Top()) {
    Printf("INCONSISTENT LOCATION PATH PAIR PRODUCED FOR TRACK %d", itr);
#ifdef VECGEOM_ROOT
    Printf("REAL");
    a->GetCurrentNode()->Print();
    Printf("REPORTED");
    fPathV[itr]->GetCurrentNode()->Print();
    //  print_trace();
#endif
  }

  // release object
  VECGEOM_NAMESPACE::NavigationState::ReleaseInstance(a);
}
#endif

#ifdef USE_VECGEOM_NAVIGATOR
GEANT_CUDA_BOTH_CODE
void GeantTrack_v::NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep,
                                              const Double_t *x, const Double_t *y,
                                              const Double_t *z, const Double_t *dirx,
                                              const Double_t *diry, const Double_t *dirz,
                                              VolumePath_t **pathin, VolumePath_t **pathout,
                                              Double_t *step, Double_t *safe, Bool_t *isonbdr,
                                              const GeantTrack_v * /*trk*/) {
  // Printf("In vec find next boundary and step\n");
  using VECGEOM_NAMESPACE::SimpleNavigator;
  using VECGEOM_NAMESPACE::Precision;
  using VECGEOM_NAMESPACE::Vector3D;
  using VECGEOM_NAMESPACE::GeoManager;
  typedef Vector3D<Precision> Vector3D_t;

  //     VolumePath_t * a = new VolumePath_t( GeoManager::Instance().getMaxDepth() );

  SimpleNavigator nav;
  for (Int_t i = 0; i < ntracks; ++i) {
    if (safe[i] > pstep[i]) {
       step[i] = pstep[i];
       isonbdr[i] = false;
       continue;
    }   
#ifdef VERBOSE
    if (pstep[i] < 0.) {
      std::cerr << " NEGATIVE PSTEP " << pstep[i] << "\n";
    }
#endif

    //    	a->Clear();
    //    	nav.LocatePoint( GeoManager::Instance().GetWorld(),
    //    			Vector3D_t( x[i], y[i], z[i] ), *a, true );
    //        if( a->Top() != NULL && a->Top() != pathin[i]->Top() )
    //         {
    //             Printf("INCONSISTENT PATH TRACK %d, boundary state %d", i, isonbdr[i] );
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
    //                                         Vector3D_t( dirx[i], diry[i], dirz[i] )  /*global
    //                                         dir*/ ,
    //                                         *a);
    //         }

    //      assert( a->Top() == pathin[i]->Top() );
    nav.FindNextBoundaryAndStep(Vector3D_t(x[i], y[i], z[i]) /* global pos */,
                                Vector3D_t(dirx[i], diry[i], dirz[i]) /* global dir */, *pathin[i],
                                *pathout[i] /* the paths */, Math::Min(1.E20, pstep[i]), step[i]);
    step[i] = Math::Max(2. * gTolerance, step[i]);
    safe[i] = (isonbdr[i]) ? 0 : nav.GetSafety(Vector3D_t(x[i], y[i], z[i]), *pathin[i]);
    safe[i] = (safe[i] < 0) ? 0. : safe[i];

#ifdef CROSSCHECK
    //************
    // CROSS CHECK USING TGEO
    //************
    TGeoNavigator *rootnav = gGeoManager->GetCurrentNavigator();
    rootnav->ResetState();
    rootnav->SetCurrentPoint(x[i], y[i], z[i]);
    rootnav->SetCurrentDirection(dirx[i], diry[i], dirz[i]);
    TGeoBranchArray *tmp = pathin[i]->ToTGeoBranchArray();
    tmp->UpdateNavigator(rootnav);
    delete tmp;
    rootnav->FindNextBoundaryAndStep(Math::Min(1.E20, pstep[i]), !isonbdr[i]);
    double stepcmp = Math::Max(2 * gTolerance, rootnav->GetStep());
    double safecmp = rootnav->GetSafeDistance();
    // pathin[i]->GetCurrentNode()->Print();
    // Printf("## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[i], step[i], stepcmp);
    // Printf("## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf",
    // pstep[i],
    //       isonbdr[i], safe[i], rootnav->Safety());

    // check nextpath
    tmp = pathout[i]->ToTGeoBranchArray();
    tmp->InitFromNavigator(rootnav);
    // Printf("## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", pathout[i]->GetCurrentNode(),
    // tmp->GetCurrentNode());
    // Printf("## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", pathout[i]->IsOnBoundary(),
    // rootnav->IsOnBoundary());

    // if( safe[i] != safecmp )
    // {
    //    nav.InspectSafetyForPoint(
    //                               Vector3D_t( x[i], y[i], z[i] )  /*global pos*/,
    //                               *pathin[i] );
    //}
    if (std::fabs(step[i] - stepcmp) > 1E-6) {
      Printf("## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[i], step[i], stepcmp);
      Printf("## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT %lf",
             pstep[i], isonbdr[i], safe[i], rootnav->Safety());

      // check nextpath
      tmp = pathout[i]->ToTGeoBranchArray();
      tmp->InitFromNavigator(rootnav);
      Printf("## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", pathout[i]->GetCurrentNode(),
             tmp->GetCurrentNode());
      Printf("## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", pathout[i]->IsOnBoundary(),
             rootnav->IsOnBoundary());

      Printf("INCONSISTENT STEP");
      nav.InspectEnvironmentForPointAndDirection(Vector3D_t(x[i], y[i], z[i]) /*global pos*/,
                                                 Vector3D_t(dirx[i], diry[i], dirz[i]), *pathin[i]);
    }
//    if( pathout[i]->IsOnBoundary() != rootnav->IsOnBoundary() )
//          {
//           Printf("INCONSISTENT BOUNDARY");
//           Printf("## PSTEP %lf VECGEOMSTEP %lf ROOTSTEP %lf", pstep[i], step[i], stepcmp);
//                 Printf("## PSTEP %lf ONBOUND %d VECGEOMSAFETY %lf ROOTSAFETY %lf BRUTEFORCEROOT
//                 %lf", pstep[i],
//                            isonbdr[i], safe[i], rootnav->Safety());
//
//                     // check nextpath
//                     tmp = pathout[i]->ToTGeoBranchArray();
//                     tmp->InitFromNavigator( rootnav );
//                     Printf("## VECGEOMNEXTNODE %p ROOTNEXTNODE %p", pathout[i]->GetCurrentNode(),
//                     tmp->GetCurrentNode());
//                    Printf("## VECGEOMBOUNDARY %d ROOTBOUNDARY %d", pathout[i]->IsOnBoundary(),
//                    rootnav->IsOnBoundary());
//
//           // nav.InspectEnvironmentForPointAndDirection(
//           //           Vector3D_t( x[i], y[i], z[i] )  /*global pos*/ ,
//            //          Vector3D_t( dirx[i], diry[i], dirz[i] ),
//             //                             *pathin[i] );
//          }
#endif
    // onboundary with respect to new point
    isonbdr[i] = pathout[i]->IsOnBoundary();

#ifdef VERBOSE
    Printf("navfindbound on %p track %d with pstep %lf yields step %lf and safety %lf\n", this, i,
           pstep[i], step[i], safe[i]);
#endif
  }
  //   delete a;
}
#else
//______________________________________________________________________________
void GeantTrack_v::NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep,
                                              const Double_t *x, const Double_t *y,
                                              const Double_t *z, const Double_t *dirx,
                                              const Double_t *diry, const Double_t *dirz,
                                              VolumePath_t **pathin, VolumePath_t **pathout,
                                              Double_t *step, Double_t *safe, Bool_t *isonbdr,
                                              const GeantTrack_v * /*trk*/) {
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
  for (Int_t i = 0; i < ntracks; i++) {
//    if (fPstepV[i] < 1.E-10) {
      // Printf("Error pstep");
//    }
    // Check if current safety allows for pstep
    if (safe[i] > pstep[i]) {
       step[i] = pstep[i];
       isonbdr[i] = false;
       continue;
    }   
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
    //      nav->SetLastSafetyForPoint(safe[i], x[i], y[i], z[i]);
    nav->FindNextBoundaryAndStep(Math::Min(1.E20, pstep[i]), !isonbdr[i]);
    step[i] = Math::Max(2 * gTolerance, nav->GetStep());
    safe[i] = isonbdr[i] ? 0. : nav->GetSafeDistance();
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
    VECGEOM_NAMESPACE::NavigationState vecgeom_in_state(
        VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
    VECGEOM_NAMESPACE::NavigationState vecgeom_out_state(
        VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
    vecgeom_in_state = *pathin[i];
    VECGEOM_NAMESPACE::SimpleNavigator vecnav;
    double vecgeom_step;
    typedef VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision> Vector3D_t;
    vecnav.FindNextBoundaryAndStep(Vector3D_t(x[i], y[i], z[i]) /* global pos */,
                                   Vector3D_t(dirx[i], diry[i], dirz[i]) /* global dir */,
                                   vecgeom_in_state, vecgeom_out_state /* the paths */,
                                   Math::Min(1.E20, pstep[i]), vecgeom_step);
    vecgeom_step = Math::Max(2 * gTolerance, vecgeom_step);
    double vecgeom_safety;
    vecgeom_safety = vecnav.GetSafety(Vector3D_t(x[i], y[i], z[i]), vecgeom_in_state);
    vecgeom_safety = (vecgeom_safety < 0) ? 0. : vecgeom_safety;
    Printf("--VECGEOM-- TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg TOBOUND %d", i, isonbdr[i],
           pstep[i], vecgeom_step, vecgeom_safety, vecgeom_out_state.IsOnBoundary());
// end crosscheck with what VECGEOM WOULD GIVE IN THIS SITUATION
// ---------------------------------------------------------
#endif

    isonbdr[i] = nav->IsOnBoundary();
  }
}
#endif

//______________________________________________________________________________
void GeantTrack_v::NavIsSameLocation(Int_t ntracks, VolumePath_t **start, VolumePath_t **end,
                                     Bool_t *same) {
  // Implementation of TGeoNavigator::IsSameLocation with vector input
  for (Int_t i = 0; i < ntracks; i++) {
    same[i] = NavIsSameLocationSingle(i, start, end);
  }
}

#ifdef USE_VECGEOM_NAVIGATOR
//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
Bool_t GeantTrack_v::NavIsSameLocationSingle(Int_t itr, VolumePath_t **start, VolumePath_t **end) {
#ifdef VERBOSE
  Printf("In NavIsSameLocation single %p for track %d", this, itr);
#endif
  // TODO: We should provide this function as a static function
  VECGEOM_NAMESPACE::SimpleNavigator simplenav;

  // this creates a tmpstate and copies in the state from end[itr]
  // we should avoid the creation of a state object here and rather use
  // some thread data?
  // was: VECGEOM_NAMESPACE::NavigationState tmpstate( *end[itr] );
  // new:
  VECGEOM_NAMESPACE::NavigationState *tmpstate =
      VECGEOM_NAMESPACE::NavigationState::MakeInstance(end[itr]->GetMaxLevel());

// cross check with answer from ROOT
#ifdef CROSSCHECK
  TGeoBranchArray *sb = start[itr]->ToTGeoBranchArray();
  TGeoBranchArray *eb = end[itr]->ToTGeoBranchArray();
#endif

  // TODO: not using the direction yet here !!
  bool samepath = simplenav.HasSamePath(VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision>(
                                            fXposV[itr], fYposV[itr], fZposV[itr]),
                                        *start[itr], *tmpstate);
  if (!samepath) {
    // Printf("CORRECTING STATE FOR TRACK %d", itr);
    // start[itr]->GetCurrentNode()->Print();
    tmpstate->CopyTo(end[itr]);
    // end[itr]->GetCurrentNode()->Print();
    // assert(end[itr]->Top() != start[itr]->Top());
  }

#ifdef CROSSCHECK
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
  nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  sb->UpdateNavigator(nav);
  bool rootsame = nav->IsSameLocation(fXposV[itr], fYposV[itr], fZposV[itr], kTRUE);
  if (rootsame != samepath) {
    Printf("INCONSISTENT ANSWER ROOT(%d) VECGEOM(%d)", rootsame, samepath);
    std::cout << VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision>(fXposV[itr], fYposV[itr],
                                                                           fZposV[itr]) << "\n";
    Printf("old state");
    sb->Print();
    nav->ResetState();
    nav->SetLastSafetyForPoint(0, 0, 0, 0);
    nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
    nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
    sb->UpdateNavigator(nav);
    nav->InspectState();
    bool rootsame = nav->IsSameLocation(fXposV[itr], fYposV[itr], fZposV[itr], kTRUE);
    nav->InspectState();
    eb->InitFromNavigator(nav);
    Printf("new state");
    eb->Print();
    Printf("VERSUS VECGEOM OLD AND NEW");
    start[itr]->printVolumePath();
    end[itr]->printVolumePath();
  } else {
    //  Printf("CONSISTENT SAME LOCATION");
  }

  delete sb;
  delete eb;
#endif // CROSSCHECK
  VECGEOM_NAMESPACE::NavigationState::ReleaseInstance(tmpstate);

  return samepath;
}
#else
//______________________________________________________________________________
Bool_t GeantTrack_v::NavIsSameLocationSingle(Int_t itr, VolumePath_t **start, VolumePath_t **end) {
  // Implementation of TGeoNavigator::IsSameLocation for single particle
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetLastSafetyForPoint(0, 0, 0, 0);
  nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
  nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  start[itr]->UpdateNavigator(nav);
  if (!nav->IsSameLocation(fXposV[itr], fYposV[itr], fZposV[itr], kTRUE)) {
    end[itr]->InitFromNavigator(nav);
    return kFALSE;
  }
  return kTRUE;
}
#endif

//______________________________________________________________________________
Int_t GeantTrack_v::SortByStatus(TrackStatus_t status) {
  // Sort tracks by a given status.
  Int_t nsel = 0;
  Int_t ntracks = GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (fStatusV[itr] == status) {
      Select(itr);
      nsel++;
    }
  }
  if (nsel) {
    if (nsel < ntracks)
      Reshuffle();
    else
      DeselectAll();
  }
  return nsel;
}


//______________________________________________________________________________
Int_t GeantTrack_v::SortByLimitingDiscreteProcess() {
  // Sort tracks for which the step was limited by discrete processes.
  Int_t nsel = 0;
  Int_t ntracks = GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (fStatusV[itr] == kPhysics && fEindexV[itr]==1000) {
      Select(itr);
      nsel++;
    }
  }
  if (nsel) {
    if (nsel < ntracks)
      Reshuffle();
    else
      DeselectAll();
  }
  return nsel;
}


//______________________________________________________________________________
Int_t GeantTrack_v::RemoveByStatus(TrackStatus_t status, GeantTrack_v &output) {
  // Remove tracks with given status from the container to the output vector,
  // then compact.
  Int_t nremoved = 0;
  Int_t ntracks = GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (fStatusV[itr] == status) {
      MarkRemoved(itr);
      nremoved++;
    }
  }
  if (!fCompact)
    Compact(&output);
  return nremoved;
}

//______________________________________________________________________________
void GeantTrack_v::PrintTrack(Int_t itr) const {
  // Print info for a given track
  const char *status[8] = {"alive",     "killed",  "inflight",  "boundary",
                           "exitSetup", "physics", "postponed", "new"};
#ifdef USE_VECGEOM_NAVIGATOR
  printf(
      "Object %p, Track %d: evt=%d slt=%d part=%d pdg=%d g5c=%d chg=%d proc=%d vid=%d nstp=%d spc=%d status=%s mass=%g\
              xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g bdr=%d\n\n",
      (const void *)this, itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fPDGV[itr],
      fG5codeV[itr], fChargeV[itr], fProcessV[itr], fVindexV[itr], fNstepsV[itr],
      (Int_t)fSpeciesV[itr], status[Int_t(fStatusV[itr])], fMassV[itr], fXposV[itr], fYposV[itr],
      fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr], fTimeV[itr],
      fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fFrombdrV[itr]);

#else
  TString path;
  fPathV[itr]->GetPath(path);
  TString nextpath;
  fNextpathV[itr]->GetPath(nextpath);

  printf("Track %d: evt=%d slt=%d part=%d pdg=%d g5c=%d eind=%d chg=%d proc=%d vid=%d nstp=%d "
         "spc=%d status=%s mass=%g xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g "
         "time=%g edep=%g pstp=%g stp=%g snxt=%g saf=%g bdr=%d\n pth=%s npth=%s\n",
         itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fPDGV[itr], fEindexV[itr],
         fG5codeV[itr], fChargeV[itr], fProcessV[itr], fVindexV[itr], fNstepsV[itr],
         (Int_t)fSpeciesV[itr], status[Int_t(fStatusV[itr])], fMassV[itr], fXposV[itr], fYposV[itr],
         fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr], fTimeV[itr],
         fEdepV[itr], fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fFrombdrV[itr],
         path.Data(), nextpath.Data());
#endif
}

//______________________________________________________________________________
void GeantTrack_v::PrintTracks() const {
  // Print all tracks
  Int_t ntracks = GetNtracks();
  for (Int_t i = 0; i < ntracks; i++)
    PrintTrack(i);
}

#ifdef USE_VECGEOM_NAVIGATOR

GEANT_CUDA_BOTH_CODE
void GeantTrack_v::ComputeTransportLength(Int_t ntracks) {
#ifndef GEANT_CUDA_DEVICE_BUILD
  static std::atomic<int> icalls(0);
  ++icalls;
#endif
  Int_t itr;
  // call the vector interface of GeantTrack_v
  NavFindNextBoundaryAndStep(ntracks, fPstepV, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV,
                             fPathV, fNextpathV, fSnextV, fSafetyV, fFrombdrV, this);
  // now we should have updated everything

  // perform a couple of additional checks/ set status flags and so on
  for (itr = 0; itr < ntracks; ++itr) {
    if ((fNextpathV[itr]->IsOutside() && fSnextV[itr] < 1.E-6) || fSnextV[itr] > 1.E19)
      fStatusV[itr] = kExitingSetup;
  }
}
#else
//______________________________________________________________________________
void GeantTrack_v::ComputeTransportLength(Int_t ntracks) {
  // Computes snext and safety for an array of tracks. For charged tracks these are the only
  // computed values, while for neutral ones the next node is checked and the boundary flag is set
  // if
  // closer than the proposed physics step.

  for (Int_t i = 0; i < ntracks; ++i) {
    ComputeTransportLengthSingle(i);
  }
}
#endif

#ifdef USE_VECGEOM_NAVIGATOR
GEANT_CUDA_BOTH_CODE
void GeantTrack_v::ComputeTransportLengthSingle(Int_t itr) {
// Computes snext and safety for a single track. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.
#ifndef GEANT_CUDA_DEVICE_BUILD
  static std::atomic<int> icalls(0);
  ++icalls;
#endif
  // inits navigator with current state
  // TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  // nav->ResetState();
  // nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
  // nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  // fPathV[itr]->UpdateNavigator(nav);
  // nav->SetLastSafetyForPoint(fSafetyV[itr], fXposV[itr], fYposV[itr], fZposV[itr]);
  // nav->FindNextBoundaryAndStep( Math::Min(1.E20, fPstepV[itr]), !fFrombdrV[itr] );

  //
  using VECGEOM_NAMESPACE::SimpleNavigator;
  using VECGEOM_NAMESPACE::Precision;
  using VECGEOM_NAMESPACE::Vector3D;
  typedef Vector3D<Precision> Vector3D_t;
/*
  if (fPstepV[itr] < fSafetyV[itr]) {
    fSnextV[itr] = fPstepV[itr];
    *fNextpathV[itr] = *fPathV[itr];
    fFrombdrV[itr] = false;
    return;
  }*/
  VECGEOM_NAMESPACE::SimpleNavigator nav;
  double step;
  nav.FindNextBoundaryAndStep(Vector3D_t(fXposV[itr], fYposV[itr], fZposV[itr]),
                              Vector3D_t(fXdirV[itr], fYdirV[itr], fZdirV[itr]), *fPathV[itr],
                              *fNextpathV[itr], Math::Min(1.E20, fPstepV[itr]), step);

  // get back step, safety, new geometry path, and other navigation information
  fSnextV[itr] = Math::Max(2 * gTolerance, step);
  fSafetyV[itr] =
      (fFrombdrV[itr]) ? 0 : nav.GetSafety(Vector3D_t(fXposV[itr], fYposV[itr], fZposV[itr]),
                                           *fPathV[itr]);
  fSafetyV[itr] = (fSafetyV[itr] < 0) ? 0. : fSafetyV[itr];
  fFrombdrV[itr] = fNextpathV[itr]->IsOnBoundary();

  // if outside detector or enormous step mark particle as exiting the detector
  if (fNextpathV[itr]->IsOutside() || fSnextV[itr] > 1.E19)
    fStatusV[itr] = kExitingSetup;

  // force track to cross under certain conditions
}
#else
void GeantTrack_v::ComputeTransportLengthSingle(Int_t itr) {
  // Computes snext and safety for a single track. For charged tracks these are the only
  // computed values, while for neutral ones the next node is checked and the boundary flag is set
  // if closer than the proposed physics step.

  // In case the proposed step is within safety, no need to compute distance to next boundary
  if (fPstepV[itr] < fSafetyV[itr]) {
    fSnextV[itr] = fPstepV[itr];
    *fNextpathV[itr] = *fPathV[itr];
    fFrombdrV[itr] = false;
    return;
  }
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  nav->ResetState();
  nav->SetCurrentPoint(fXposV[itr], fYposV[itr], fZposV[itr]);
  nav->SetCurrentDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  fPathV[itr]->UpdateNavigator(nav);
  //   nav->SetLastSafetyForPoint(fSafetyV[itr], fXposV[itr], fYposV[itr], fZposV[itr]);
  nav->FindNextBoundaryAndStep(Math::Min(1.E20, fPstepV[itr]), !fFrombdrV[itr]);
  fSnextV[itr] = Math::Max(2 * gTolerance, nav->GetStep());
  fSafetyV[itr] = nav->GetSafeDistance();
  fNextpathV[itr]->InitFromNavigator(nav);

#ifdef CROSSCHECK
  VECGEOM_NAMESPACE::NavigationState vecgeom_in_state(
      VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
  VECGEOM_NAMESPACE::NavigationState vecgeom_out_state(
      VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
  vecgeom_in_state = *fPathV[itr];
  VECGEOM_NAMESPACE::SimpleNavigator vecnav;
  double vecgeom_step;
  typedef VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision> Vector3D_t;
  vecnav.FindNextBoundaryAndStep(Vector3D_t(fXposV[itr], fYposV[itr], fZposV[itr]) /* global pos */,
                                 Vector3D_t(fXdirV[itr], fYdirV[itr], fZdirV[itr]) /* global dir */,
                                 vecgeom_in_state, vecgeom_out_state /* the paths */,
                                 Math::Min(1.E20, fPstepV[itr]), vecgeom_step);
  vecgeom_step = Math::Max(2 * gTolerance, vecgeom_step);
  double vecgeom_safety;
  vecgeom_safety =
      vecnav.GetSafety(Vector3D_t(fXposV[itr], fYposV[itr], fZposV[itr]), vecgeom_in_state);
  vecgeom_safety = (vecgeom_safety < 0) ? 0. : vecgeom_safety;
  //   Printf("--VECGEOM-- TRACK %d BOUND %d PSTEP %lg STEP %lg SAFETY %lg TOBOUND %d",
  //                  i, isonbdr[i], pstep[i], vecgeom_step, vecgeom_safety,
  //                  vecgeom_out_state.IsOnBoundary());
  if (!(vecgeom_out_state.IsOutside() || fNextpathV[itr]->IsOutside()) &&
      (vecgeom_out_state.GetCurrentNode() != fNextpathV[itr]->GetCurrentNode())) {
    warnings++;
    InspectGeometryState(itr);
    Warning("", "NavigationWarning");
    if (warnings > 100)
      assert(!"NOT SAME RESULT FOR NEXTNODE");
  }
#endif
  fFrombdrV[itr] = nav->IsOnBoundary();

  // if outside detector or enormous step mark particle as exiting the detector
  if (fNextpathV[itr]->IsOutside() || fSnextV[itr] > 1.E19)
    fStatusV[itr] = kExitingSetup;
}
#endif
//______________________________________________________________________________
TransportAction_t GeantTrack_v::PostponedAction(Int_t ntracks) const {
  // Check the action to be taken according the current policy
  if (!ntracks)
    return kDone;
  // Temporary hook
  if (ntracks == 1) {
    //      if (gPropagator->GetPolicy()<GeantPropagator::kPrioritize)
    //           return kPostpone;
    //      else
    return kSingle;
  }
  return kVector;
}

//______________________________________________________________________________
Int_t GeantTrack_v::PropagateTracks(GeantTrack_v &output, Int_t tid) {
  // Propagate the ntracks in the current volume with their physics steps (already
  // computed)
  // Vectors are pushed downstream when efficient.
  Int_t ntracks = GetNtracks();
  // Check if tracking the remaining tracks can be postponed
  TransportAction_t action = PostponedAction(ntracks);
  if (action == kPostpone) {
    PostponeTracks(output);
    return 0;
  }
  if (action != kVector)
    return PropagateTracksScalar(output, tid, 0);
  // Compute transport length in geometry, limited by the physics step
  ComputeTransportLength(ntracks);
  //         Printf("====== After ComputeTransportLength:");
  //         PrintTracks();

  Int_t itr = 0;
  Int_t icrossed = 0;
  Int_t nsel = 0;
  Double_t lmax;
  const Double_t eps = 1.E-2; // 100 micron
  const Double_t bmag = gPropagator->fBmag;
  //   TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  //   Int_t tid = nav->GetThreadId();
  //   GeantThreadData *td = gPropagator->fThreadData[tid];
  // Remove dead tracks, propagate neutrals
  for (itr = 0; itr < ntracks; itr++) {
    // Mark dead tracks for copy/removal
    if (fStatusV[itr] == kKilled) {
      MarkRemoved(itr);
      continue;
    }
    // Propagate straight tracks to the precomputed location and update state,
    // then mark them for copy/removal
    // (Inlined from PropagateStraight)
    if (fChargeV[itr] == 0 || bmag < 1.E-10) {
      // Do straight propagation to physics process or boundary
      if (fFrombdrV[itr]) {
        // *fPathV[itr] = *fNextpathV[itr];
        if (fNextpathV[itr]->IsOutside())
          fStatusV[itr] = kExitingSetup;
        else
          fStatusV[itr] = kBoundary;
        icrossed++;
      } else {
        fStatusV[itr] = kPhysics;
      }
      fPstepV[itr] -= fSnextV[itr];
      fStepV[itr] += fSnextV[itr];
      fSafetyV[itr] -= fSnextV[itr];
      if (fSafetyV[itr] < 0.)
        fSafetyV[itr] = 0;
      fXposV[itr] += fSnextV[itr] * fXdirV[itr];
      fYposV[itr] += fSnextV[itr] * fYdirV[itr];
      fZposV[itr] += fSnextV[itr] * fZdirV[itr];
      fSnextV[itr] = 0;
      fNstepsV[itr]++;
#ifndef GEANT_CUDA_DEVICE_BUILD
      gPropagator->fNsnextSteps++; // should use atomics
#endif
      MarkRemoved(itr);
#ifdef USE_VECGEOM_NAVIGATOR
//            CheckLocationPathConsistency(itr);
#endif
    }
  }
  // Compact remaining tracks and move the removed oned to the output container
  if (!fCompact)
    Compact(&output);
  // Check if tracking the remaining tracks can be postponed
  action = PostponedAction(fNtracks);
  switch (action) {
  case kDone:
    return icrossed;
  case kSingle:
    icrossed += PropagateTracksScalar(output, tid, 1);
    return icrossed;
  case kPostpone:
    PostponeTracks(output);
    return icrossed;
  case kVector:
    break;
  }
  // REMAINING ONLY CHARGED TRACKS IN MAG. FIELD
  // Continue with vectorized mode ...

  // New algorithm: we use the track sagitta to estimate the "bending" error,
  // i.e. what is the propagated length for which the track deviation in magnetic
  // field with respect to straight propagation is less than epsilon.
  // Take the maximum between the safety and the "bending" safety
  nsel = 0;
  ntracks = GetNtracks();
  Double_t *steps = GeantPropagator::Instance()->fThreadData[tid]->GetDblArray(ntracks);
  for (itr = 0; itr < fNtracks; itr++) {
    lmax = SafeLength(itr, eps);
    lmax = Math::Max(lmax, fSafetyV[itr]);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if frombdr=1) or the proposed  physics
    // step (frombdr=0)
    steps[itr] = (fFrombdrV[itr]) ? Math::Min(lmax, TMath::Max(fSnextV[itr],1.E-4))
                                  : Math::Min(lmax, fPstepV[itr]);
//    if (fFrombdrV[itr] && steps[itr]<1.E-8) steps[itr] = 1.E-3;
    //Printf("snext=%g lmax=%g", fSnextV[itr], lmax);
    //      Printf("track %d: step=%g (safelen=%g)", itr, steps[itr], lmax);
  }
  // Propagate the vector of tracks
  PropagateInVolume(ntracks, steps, tid);
  //         Printf("====== After PropagateInVolume:");
  //         PrintTracks();
  // Some tracks made it to physics steps (kPhysics)
  //         -> remove and copy to output
  // Some tracks were propagated with steps greater than safety
  //         -> check possible crossing via NavIsSameLocation
  // Some tracks were propagated with steps less than safety
  //         -> keep in the container

  // Select tracks that made it to physics and copy to output
  ntracks = GetNtracks();
  for (itr = 0; itr < ntracks; itr++) {
    if (fStatusV[itr] == kPhysics) {
      MarkRemoved(itr);
      fNstepsV[itr]++;
    }
  }
  if (!fCompact)
    Compact(&output);
  // Select tracks that are in flight or were propagated to boundary with
  // steps bigger than safety
  nsel = 0;
  ntracks = GetNtracks();
  for (itr = 0; itr < ntracks; itr++) {
    if (fSafetyV[itr] < 1.E-10 || fSnextV[itr] < 1.E-10) {
      Select(itr);
      nsel++;
    }
  }
  // The selected tracks may have crossed the boundaries - check that
  if (nsel) {
    if (nsel < GetNtracks())
      Reshuffle();
    else
      DeselectAll();
    Bool_t *same = GeantPropagator::Instance()->fThreadData[tid]->GetBoolArray(nsel);
    NavIsSameLocation(nsel, fPathV, fNextpathV, same);
    for (itr = 0; itr < nsel; itr++) {
      if (same[itr])
        continue;
      // Boundary crossed -> update current path
      //*fPathV[itr] = *fNextpathV[itr];
      fStatusV[itr] = kBoundary;
      if (fNextpathV[itr]->IsOutside())
        fStatusV[itr] = kExitingSetup;
      fFrombdrV[itr] = kTRUE;
      icrossed++;
      fNstepsV[itr]++;
      MarkRemoved(itr);
    }
    //         Printf("====== After finding crossing tracks (ncross=%d):", icrossed);
    //         PrintTracks();
    if (!fCompact)
      Compact(&output);
  }
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Int_t GeantTrack_v::PropagateSingleTrack(GeantTrack_v & /*output*/, Int_t itr, Int_t tid, Int_t stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  Int_t icrossed = 0;
  Double_t step, lmax;
  const Double_t eps = 1.E-4; // 1 micron
#ifdef GEANT_CUDA_DEVICE_BUILD
  const Double_t bmag = gPropagator_fBmag;
#else
  const Double_t bmag = gPropagator->fBmag;
#endif

  // Mark dead tracks for copy/removal
  if (fStatusV[itr] == kKilled) {
     MarkRemoved(itr);
     return icrossed;
  }
  // Compute transport length in geometry, limited by the physics step
  ComputeTransportLengthSingle(itr);
  //      Printf("====== After ComputeTransportLengthSingle:");
  //      PrintTrack(itr);
  // Stage 0: straight propagation
  if (stage == 0) {
     if (fChargeV[itr] == 0 || bmag < 1.E-10) {
        // Do straight propagation to physics process or boundary
        if (fFrombdrV[itr]) {
          //*fPathV[itr] = *fNextpathV[itr];
          if (fNextpathV[itr]->IsOutside())
            fStatusV[itr] = kExitingSetup;
          else
            fStatusV[itr] = kBoundary;
          icrossed++;
        } else {
           fStatusV[itr] = kPhysics;
        }
        fPstepV[itr] -= fSnextV[itr];
        fStepV[itr] += fSnextV[itr];
        fSafetyV[itr] -= fSnextV[itr];
        if (fSafetyV[itr] < 0.)
           fSafetyV[itr] = 0;
        fXposV[itr] += fSnextV[itr] * fXdirV[itr];
        fYposV[itr] += fSnextV[itr] * fYdirV[itr];
        fZposV[itr] += fSnextV[itr] * fZdirV[itr];
        fSnextV[itr] = 0;
        fNstepsV[itr]++;
#ifndef GEANT_CUDA_DEVICE_BUILD
        gPropagator->fNsnextSteps++;
#endif
        MarkRemoved(itr);
#ifdef USE_VECGEOM_NAVIGATOR
        //            CheckLocationPathConsistency(itr);
#endif
        return icrossed;
     }
  }
  // Stage 1: mag field propagation for tracks with pstep<safety
  if (stage <= 1) {
     // REMAINING ONLY CHARGED TRACKS IN MAG. FIELD
     // New algorithm: we use the track sagitta to estimate the "bending" error,
     // i.e. what is the propagated length for which the track deviation in magnetic
     // field with respect to straight propagation is less than epsilon.
     // Take the maximum between the safety and the "bending" safety
     lmax = SafeLength(itr, eps);
     lmax = Math::Max(lmax, fSafetyV[itr]);
     // Select step to propagate as the minimum among the "safe" step and:
     // the straight distance to boundary (if frombdr=1) or the proposed  physics
     // step (frombdr=0)
     step = (fFrombdrV[itr]) ? Math::Min(lmax, fSnextV[itr] + 10 * gTolerance)
        : Math::Min(lmax, fPstepV[itr]);
     //      Printf("track %d: step=%g (safelen=%g)", itr, step, lmax);
     PropagateInVolumeSingle(itr, step, tid);
     //      Printf("====== After PropagateInVolumeSingle:");
     //      PrintTrack(itr);
     // The track may have made it to physics steps (kPhysics)
     //         -> remove and copy to output
     // The track may have been propagated with step greater than safety
     //         -> check possible crossing via NavIsSameLocation
     // The track may have been propagated with step less than safety
     //         -> keep in the container

     // Select tracks that made it to physics and copy to output
     if (fStatusV[itr] == kPhysics) {
        MarkRemoved(itr);
        fNstepsV[itr]++;
        return icrossed;
     }
     // Select tracks that are in flight or were propagated to boundary with
     // steps bigger than safety
     if (fSafetyV[itr] < 1.E-10 || fSnextV[itr] < 1.E-10) {
        Bool_t same = NavIsSameLocationSingle(itr, fPathV, fNextpathV);
        if (same)
           return icrossed;
        // Boundary crossed -> update current path
        //*fPathV[itr] = *fNextpathV[itr];
        fStatusV[itr] = kBoundary;
        if (fNextpathV[itr]->IsOutside())
          fStatusV[itr] = kExitingSetup;
        fFrombdrV[itr] = kTRUE;
        icrossed++;
        fNstepsV[itr]++;
        MarkRemoved(itr);
     }
#ifdef USE_VECGEOM_NAVIGATOR
     //         CheckLocationPathConsistency(itr);
#endif
  }
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Int_t GeantTrack_v::PropagateTracksScalar(GeantTrack_v &output, Int_t tid, Int_t stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  Int_t icrossed = 0;
  Int_t ntracks = GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
     icrossed += PropagateSingleTrack(output, itr, tid, stage);
  }
//   Printf("====== After finding crossing tracks (ncross=%d):", icrossed);
//   PrintTracks();
// Compact remaining tracks and move the removed oned to the output container
#ifndef GEANT_CUDA_DEVICE_BUILD
  if (!fCompact)
    Compact(&output);
#endif
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Double_t GeantTrack_v::Curvature(Int_t i) const {
  // Curvature assuming constant field is along Z
  const Double_t tiny = 1.E-30;
#ifdef GEANT_CUDA_DEVICE_BUILD
  const Double_t bmag = gPropagator_fBmag;
#else
  const Double_t bmag = gPropagator->fBmag;
#endif
  return Math::Abs(kB2C * fChargeV[i] * bmag / (Pt(i) + tiny));
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Double_t GeantTrack_v::SafeLength(Int_t i, Double_t eps) {
  // Returns the propagation length in field such that the propagated point is
  // shifted less than eps with respect to the linear propagation.
  Double_t c = Curvature(i);
  if (c < 1.E-10)
    return 1.E20;
  return 2. * Math::Sqrt(eps / c);
}

//______________________________________________________________________________
Int_t GeantTrack_v::PostponeTracks(GeantTrack_v &output) {
  // Postpone transport of remaining tracks and copy them to the output.
  Int_t npostponed = GetNtracks();
  for (Int_t itr = 0; itr < npostponed; itr++)
    fStatusV[itr] = kPostponed;
  // Move these tracks to the output container
  output.AddTracks(*this, 0, npostponed - 1, kTRUE);
  RemoveTracks(0, npostponed - 1);
  Clear();
  return npostponed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
Int_t GeantTrack_v::PostponeTrack(Int_t itr, GeantTrack_v &output) {
  // Postpone transport of a track and copy it to the output.
  // Returns where in the output the track was added.

  fStatusV[itr] = kPostponed;
  // Move these tracks to the output container
  Int_t new_itr = output.AddTrack(*this, itr, kTRUE);
  MarkRemoved(itr);
  return new_itr;
}

//______________________________________________________________________________
TGeoVolume *GeantTrack_v::GetVolume(Int_t i) const {
  // Current volume the track is into
  return ((TGeoVolume*)gGeoManager->GetListOfVolumes()->At(fVindexV[i]));
}

//______________________________________________________________________________
TGeoVolume *GeantTrack_v::GetNextVolume(Int_t i) const {
  // Next volume the track is getting into
  return fNextpathV[i]->GetCurrentNode()->GetVolume();
}

//______________________________________________________________________________
TGeoMaterial *GeantTrack_v::GetMaterial(Int_t i) const {
  // Current material the track is into
  TGeoMedium *med = GetVolume(i)->GetMedium();
  if (!med)
    return 0;
  return med->GetMaterial();
}
