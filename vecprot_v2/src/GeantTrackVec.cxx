#include "GeantTrack.h"

#include "globals.h"
#include "Geant/Error.h"
#include <execinfo.h>

#ifdef USE_VECGEOM_NAVIGATOR
#include "ScalarNavInterfaceVG.h"
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/ABBoxNavigator.h"
#include "volumes/PlacedVolume.h" // equivalent of TGeoNode
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"
#include "base/SOA3D.h"
#ifdef CROSSCHECK
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif
#else
#include "ScalarNavInterfaceTGeo.h"
#include <iostream>
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif

#include "WorkloadManager.h"

#include "GeantTaskData.h"
#include "StepChecker.h"
#include "ConstBzFieldHelixStepper.h"
#include "ConstFieldHelixStepper.h"
#include "GeantScheduler.h"

// #include "VVectorField.h"
#include "Units.h"     //  Field Units - to be 'unified'
#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "FieldLookup.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h>
#else
#include "mm_malloc.h"
#endif
#include <cassert>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

namespace host_constant {
#ifdef USE_VECGEOM_NAVIGATOR
const double gTolerance = vecgeom::kTolerance;
#else
const double gTolerance = TGeoShape::Tolerance();
#endif
}
#ifdef CUDA_SEP_COMP
namespace device_constant {
   __constant__ double gTolerance;
}
#endif

#ifdef USE_VECGEOM_NAVIGATOR
using namespace VECGEOM_NAMESPACE;
#endif

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v()
    : fNtracks(0), fNselected(0), fCompact(true), fMixed(false), fMaxtracks(0), fHoles(0), fSelected(0), fMaxDepth(0),
      fBufSize(0), fVPstart(0), fBuf(0), fEventV(0), fEvslotV(0), fParticleV(0), fMotherV(0), fPDGV(0), fGVcodeV(0), fEindexV(0), fBindexV(0),
      fChargeV(0), fProcessV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0),
      fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0), fStepV(0),
      fSnextV(0), fSafetyV(0), fNintLenV(0), fIntLenV(0), fBoundaryV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
  // Dummy ctor.
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(int size, int maxdepth)
    : fNtracks(0), fNselected(0), fCompact(true), fMixed(false), fMaxtracks(0), fHoles(0), fSelected(0),
      fMaxDepth(maxdepth), fBufSize(0), fVPstart(0), fBuf(0), fEventV(0), fEvslotV(0), fParticleV(0), fPDGV(0),
      fGVcodeV(0), fEindexV(0), fChargeV(0), fProcessV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0),
      fMassV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fNintLenV(0), fIntLenV(0), fBoundaryV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
  // Constructor with maximum capacity.
  Resize(size);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantTrack_v *GeantTrack_v::MakeInstanceAt(void *addr, unsigned int nTracks, int maxdepth) {
  return new (addr) GeantTrack_v(addr, nTracks, maxdepth);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantTrack_v::GeantTrack_v(void *addr, unsigned int nTracks, int maxdepth)
    : fNtracks(0), fNselected(0), fCompact(true), fMixed(false), fMaxtracks(GeantTrack::round_up_align(nTracks)), fHoles(0),
      fSelected(0), fMaxDepth(maxdepth), fBufSize(0), fVPstart(0), fBuf(0), fEventV(0), fEvslotV(0), fParticleV(0), fMotherV(0),
      fPDGV(0), fGVcodeV(0), fEindexV(0), fBindexV(0), fChargeV(0), fProcessV(0), fNstepsV(0), fSpeciesV(0),
      fStatusV(0), fMassV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0),
      fTimeV(0), fEdepV(0), fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fNintLenV(0), fIntLenV(0), fBoundaryV(0), fPendingV(0),
      fPathV(0), fNextpathV(0) {

  // Constructor with maximum capacity.
  fBuf = ((char *)addr) + GeantTrack::round_up_align(sizeof(GeantTrack_v));
  fBufSize = BufferSize(nTracks, maxdepth);
  memset(fBuf, 0, fBufSize);
  AssignInBuffer(fBuf, nTracks);
  memset(fPathV, 0, nTracks * sizeof(VolumePath_t *));
  memset(fNextpathV, 0, nTracks * sizeof(VolumePath_t *));
}

//______________________________________________________________________________
GeantTrack_v::GeantTrack_v(const GeantTrack_v &track_v)
    : fNtracks(0), fNselected(track_v.fNselected), fCompact(track_v.fCompact), fMixed(track_v.fMixed),
      fMaxtracks(track_v.fMaxtracks), fHoles(0), fSelected(0), fMaxDepth(track_v.fMaxDepth), fBufSize(track_v.fBufSize),
      fVPstart(0), fBuf(0), fEventV(0), fEvslotV(0), fParticleV(0), fMotherV(0), fPDGV(0), fGVcodeV(0), fEindexV(0), fBindexV(0), fChargeV(0),
      fProcessV(0), fNstepsV(0), fSpeciesV(0), fStatusV(0), fMassV(0), fXposV(0), fYposV(0), fZposV(0),
      fXdirV(0), fYdirV(0), fZdirV(0), fPV(0), fEV(0), fTimeV(0), fEdepV(0), fPstepV(0), fStepV(0), fSnextV(0),
      fSafetyV(0), fNintLenV(0), fIntLenV(0), fBoundaryV(0), fPendingV(0), fPathV(0), fNextpathV(0) {
// Copy constructor
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  fNtracks.store(track_v.fNtracks);
#else
  fNtracks = track_v.fNtracks;
#endif
  fBuf = (char *)_mm_malloc(fBufSize, GEANT_ALIGN_PADDING);
  memcpy(fBuf, track_v.fBuf, fBufSize);
  AssignInBuffer(&fBuf[0], fMaxtracks);
}

//______________________________________________________________________________
GeantTrack_v &GeantTrack_v::operator=(const GeantTrack_v &track_v) {
  // Assignment operator
  if (&track_v != this) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    fNtracks.store(track_v.fNtracks);
#else
    fNtracks = track_v.fNtracks;
#endif
    int size = track_v.fMaxtracks;
    fMaxDepth = track_v.fMaxDepth;
    fBufSize = track_v.fBufSize;
    if (fMaxtracks < size) {
      _mm_free(fBuf);
      fBuf = (char *)_mm_malloc(fBufSize, GEANT_ALIGN_PADDING);
    }
    fMaxtracks = size;
    fNselected = track_v.fNselected;
    fHoles = 0;
    fSelected = 0;
    fCompact = track_v.fCompact;
    fMixed = track_v.fMixed;
    memcpy(fBuf, track_v.fBuf, size * sizeof(GeantTrack));
    AssignInBuffer(&fBuf[0], size);
  }
  return *this;
}

//______________________________________________________________________________
GeantTrack_v::~GeantTrack_v() {
  // Destructor.
  int ntracks = GetNtracks();
  for (auto i = 0; i < ntracks; ++i) {
    VolumePath_t::ReleaseInstance(fPathV[i]);
    VolumePath_t::ReleaseInstance(fNextpathV[i]);
  }
  _mm_free(fBuf);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::AssignInBuffer(char *buff, int size) {
  // Assign all internal class arrays in the supplied buffer, padded by supplied
  // size.

  const int size_intn = size * sizeof(int);
  const int size_doublen = size * sizeof(double);
  const int size_booln = size * sizeof(bool);
  char *buf = buff;
  fEventV = (int *)buf;
  buf += size_intn;
  fEvslotV = (int *)buf;
  buf += size_intn;
  fParticleV = (int *)buf;
  buf += size_intn;
  fMotherV = (int *)buf;
  buf += size_intn;
  fPDGV = (int *)buf;
  buf += size_intn;
  fGVcodeV = (int *)buf;
  buf += size_intn;
  fEindexV = (int *)buf;
  buf += size_intn;
  fBindexV = (int *)buf;
  buf += size_intn;
  fChargeV = (int *)buf;
  buf += size_intn;
  fProcessV = (int *)buf;
  buf += size_intn;
  fNstepsV = (int *)buf;
  buf += size_intn;
  fSpeciesV = (Species_t *)buf;
  buf += size * sizeof(Species_t);
  fStatusV = (TrackStatus_t *)buf;
  buf += size * sizeof(TrackStatus_t);
  fMassV = (double *)buf;
  buf += size_doublen;
  fXposV = (double *)buf;
  buf += size_doublen;
  fYposV = (double *)buf;
  buf += size_doublen;
  fZposV = (double *)buf;
  buf += size_doublen;
  fXdirV = (double *)buf;
  buf += size_doublen;
  fYdirV = (double *)buf;
  buf += size_doublen;
  fZdirV = (double *)buf;
  buf += size_doublen;
  fPV = (double *)buf;
  buf += size_doublen;
  fEV = (double *)buf;
  buf += size_doublen;
  fTimeV = (double *)buf;
  buf += size_doublen;
  fEdepV = (double *)buf;
  buf += size_doublen;
  fPstepV = (double *)buf;
  buf += size_doublen;
  fStepV = (double *)buf;
  buf += size_doublen;
  fSnextV = (double *)buf;
  buf += size_doublen;
  fSafetyV = (double *)buf;
  buf += size_doublen;
  fNintLenV = (double *)buf;
  buf += size_doublen;
  fIntLenV = (double *)buf;
  buf += size_doublen;
  fBoundaryV = (bool *)buf;
  buf += size_booln;
  fPendingV = (bool *)buf;
  buf += size_booln;
  fPathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  fNextpathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
  fVPstart = buf;
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  // Allocate VolumePath_t objects in the reserved buffer space
  for (auto i = 0; i < 2 * size; ++i)
    VolumePath_t::MakeInstanceAt(fMaxDepth, buf + i * size_vpath);
  buf += 2 * size * size_vpath;

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
  size_t size_bits = BitSet::SizeOfInstance(size);
  fHoles = BitSet::MakeInstanceAt(size, buf);
  buf += size_bits;

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
  fSelected = BitSet::MakeInstanceAt(size, buf);
}

//______________________________________________________________________________
void GeantTrack_v::CopyToBuffer(char *buff, int size) {
  // Copy existing track arrays into new buffer, padded by supplied size
  int ntracks = GetNtracks();
  const int size_int = ntracks * sizeof(int);
  const int size_double = ntracks * sizeof(double);
  const int size_intn = size * sizeof(int);
  const int size_doublen = size * sizeof(double);
  char *buf = buff;
  memcpy(buf, fEventV, size_int);
  fEventV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fEvslotV, size_int);
  fEvslotV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fParticleV, size_int);
  fParticleV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fMotherV, size_int);
  fMotherV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fPDGV, size_int);
  fPDGV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fGVcodeV, size_int);
  fGVcodeV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fEindexV, size_int);
  fEindexV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fBindexV, size_int);
  fBindexV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fChargeV, size_int);
  fChargeV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fProcessV, size_int);
  fProcessV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fNstepsV, size_int);
  fNstepsV = (int *)buf;
  buf += size_intn;
  memcpy(buf, fSpeciesV, ntracks * sizeof(Species_t));
  fSpeciesV = (Species_t *)buf;
  buf += size * sizeof(Species_t);
  memcpy(buf, fStatusV, ntracks * sizeof(TrackStatus_t));
  fStatusV = (TrackStatus_t *)buf;
  buf += size * sizeof(TrackStatus_t);
  memcpy(buf, fMassV, size_double);
  fMassV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fXposV, size_double);
  fXposV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fYposV, size_double);
  fYposV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fZposV, size_double);
  fZposV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fXdirV, size_double);
  fXdirV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fYdirV, size_double);
  fYdirV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fZdirV, size_double);
  fZdirV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fPV, size_double);
  fPV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fEV, size_double);
  fEV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fTimeV, size_double);
  fTimeV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fEdepV, size_double);
  fEdepV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fPstepV, size_double);
  fPstepV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fStepV, size_double);
  fStepV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fSnextV, size_double);
  fSnextV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fSafetyV, size_double);
  fSafetyV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fNintLenV, size_double);
  fNintLenV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fIntLenV, size_double);
  fIntLenV = (double *)buf;
  buf += size_doublen;
  memcpy(buf, fBoundaryV, ntracks * sizeof(bool));
  fBoundaryV = (bool *)buf;
  buf += size * sizeof(bool);
  memcpy(buf, fPendingV, ntracks * sizeof(bool));
  fPendingV = (bool *)buf;
  buf += size * sizeof(bool);
  //   memcpy(buf, fPathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **pathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  //   memcpy(buf, fNextpathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **nextpathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
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

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
  size_t size_bits = BitSet::SizeOfInstance(size);
//  BitSet *holes = BitSet::MakeCopyAt(*fHoles, buf, size);
  BitSet *holes = BitSet::MakeInstanceAt(size, buf);
  BitSet::ReleaseInstance(fHoles);
  fHoles = holes;
  buf += size_bits;

  // Now the start of objects, we need to align the memory.
  buf = GeantTrack::round_up_align(buf);
  BitSet *selected = BitSet::MakeInstanceAt(size, buf);
  BitSet::ReleaseInstance(fSelected);
  fSelected = selected;
}

//______________________________________________________________________________
bool GeantTrack_v::IsSame(const GeantTrack_v &tr1, int i1, const GeantTrack_v &tr2, int i2) {
  // Compare two tracks.
  long long int chk1, chk2;
  chk1 = tr1.fEventV[i1] + tr1.fEvslotV[i1] + tr1.fParticleV[i1] + tr1.fMotherV[i1] + tr1.fPDGV[i1] + tr1.fGVcodeV[i1] + tr1.fEindexV[i1] +
         tr1.fChargeV[i1] + tr1.fProcessV[i1] + tr1.fNstepsV[i1] + (long long int)tr1.fSpeciesV[i1] +
         (long long int)tr1.fStatusV[i1];
  chk2 = tr2.fEventV[i2] + tr2.fEvslotV[i2] + tr2.fParticleV[i2] + tr2.fMotherV[i2] + tr2.fPDGV[i2] + tr2.fGVcodeV[i2] + tr2.fEindexV[i2] +
         tr2.fChargeV[i2] + tr2.fProcessV[i2] + tr2.fNstepsV[i2] + (long long int)tr2.fSpeciesV[i2] +
         (long long int)tr2.fStatusV[i2];
  if (chk1 != chk2)
    return false;
  double dchk1, dchk2;
  dchk1 = (long long int)tr1.fMassV[i1] + tr1.fXposV[i1] + tr1.fYposV[i1] + tr1.fZposV[i1] + tr1.fXdirV[i1] +
          tr1.fYdirV[i1] + tr1.fZdirV[i1] + tr1.fPV[i1] + tr1.fEdepV[i1] + tr1.fEV[i1] + tr1.fPstepV[i1] +
          tr1.fStepV[i1] + tr1.fSnextV[i1] + tr1.fSafetyV[i1] + tr1.fNintLenV[i1] + tr1.fIntLenV[i1];
  dchk2 = (long long int)tr2.fMassV[i2] + tr2.fXposV[i2] + tr2.fYposV[i2] + tr2.fZposV[i2] + tr2.fXdirV[i2] +
          tr2.fYdirV[i2] + tr2.fZdirV[i2] + tr2.fPV[i2] + tr2.fEdepV[i2] + tr2.fEV[i2] + tr2.fPstepV[i2] +
          tr2.fStepV[i2] + tr2.fSnextV[i2] + tr2.fSafetyV[i2] + tr2.fNintLenV[i2]+ tr2.fIntLenV[i2];
  if (!Math::AreEqualAbs(dchk1, dchk2, 1.E-10))
    return false;
  if (tr1.fPendingV[i1] != tr2.fPendingV[i2])
    return false;
  return true;
}

//______________________________________________________________________________
bool GeantTrack_v::IsNormalized(int itr, double tolerance) const {
  // Check if track direction is normalized within tolerance
  double norm = fXdirV[itr] * fXdirV[itr] + fYdirV[itr] * fYdirV[itr] + fZdirV[itr] * fZdirV[itr];
  if (fabs(1. - norm) > tolerance)
    return false;
  return true;
}

//______________________________________________________________________________
void GeantTrack_v::CheckTracks() {
  //  for (int i=0; i<fMaxtracks; ++i)
  //     if (fNextpathV[i]) fNextpathV[i]->SetClient(gPropagator);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
size_t GeantTrack_v::BufferSize(size_t nTracks, size_t maxdepth) {
  // return the contiguous memory size needed to hold a GeantTrack_v's data

  size_t size = GeantTrack::round_up_align(nTracks); // When called internally this ought to be a nop
  size_t size_nav = 2 * size * VolumePath_t::SizeOfInstance(maxdepth);
  // NOTE: Most likely the 'real' problem here is that BitSet::SizeOfInstance return
  // a number that is as small as possible rather than a number that is usuable to
  // be able to make array of BitSet.
  size_t size_bits = 2 * GeantTrack::round_up_align(BitSet::SizeOfInstance(size));

  // Since we already round nTracks, we only need to round the last two
  // to have the proper space for object alignment
  return size * sizeof(GeantTrack) + GeantTrack::round_up_align(size_nav) + size_bits;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
size_t GeantTrack_v::SizeOfInstance(size_t nTracks, size_t maxdepth) {
  // return the contiguous memory size needed to hold a GeantTrack_v

  return GeantTrack::round_up_align(sizeof(GeantTrack_v))+BufferSize(nTracks,maxdepth);
}

//______________________________________________________________________________
void GeantTrack_v::Resize(int newsize) {
  // Resize the container.
  int size = GeantTrack::round_up_align(newsize);
  if (size < GetNtracks()) {
    Geant::Error("Resize","%s","Cannot resize to less than current track content");
    // Geant::Error("Resize","Cannot resize to less than current track content (round up size: %d less than current number of tracks: %d)",size,GetNtracks());
    return;
  }
  fBufSize = BufferSize(size, fMaxDepth);
  if (!fCompact)
    Compact();

  char *buf = (char *)_mm_malloc(fBufSize, GEANT_ALIGN_PADDING);
  memset(buf, 0, fBufSize);
  fMaxtracks = size;
  if (!fBuf) {
    // All arrays are contiguous in a single buffer and aligned with the
    // same padding GEANT_ALIGN_PADDING
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
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::AddTrack(GeantTrack &track, bool /*import*/) {
  // Add new track to the array. If addition is done on top of non-compact array,
  // the track will be inserted without updating the number of tracks. If track is
  // imported just copy the pointers to the navigation states and reset the sources.
  // Returns the location where the track was added.
  int itrack = GetNtracks();
  if (!fCompact)
    itrack = fHoles->FirstSetBit();
  if (itrack == fMaxtracks) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    Resize(2 * fMaxtracks);
#else
    printf("Error in GeantTrack_v::AddTrack, resizing is not supported in device code\n");
#endif
  }
  fHoles->ResetBitNumber(itrack);
  fSelected->ResetBitNumber(itrack);
  fEventV[itrack] = track.Event();
  fEvslotV[itrack] = track.EventSlot();
  fParticleV[itrack] = track.Particle();
  fMotherV[itrack] = track.Mother();
  fPDGV[itrack] = track.PDG();
  fGVcodeV[itrack] = track.GVcode();
  fEindexV[itrack] = track.EIndex();
  fBindexV[itrack] = track.BIndex();
  fChargeV[itrack] = track.Charge();
  fProcessV[itrack] = track.Process();
  fNstepsV[itrack] = track.GetNsteps();
  fSpeciesV[itrack] = track.Species();
  fStatusV[itrack] = track.Status();
  fMassV[itrack] = track.Mass();
  fXposV[itrack] = track.X();
  fYposV[itrack] = track.Y();
  fZposV[itrack] = track.Z();
  fXdirV[itrack] = track.Dx();
  fYdirV[itrack] = track.Dy();
  fZdirV[itrack] = track.Dz();
  fPV[itrack] = track.P();
  fEV[itrack] = track.E();
  fTimeV[itrack] = track.Time();
  fEdepV[itrack] = track.Edep();
  fPstepV[itrack] = track.GetPstep();
  fStepV[itrack] = track.GetStep();
  fSnextV[itrack] = track.GetSnext();
  fSafetyV[itrack] = track.GetSafety();
  fNintLenV[itrack] = track.GetNintLen();
  fIntLenV[itrack] = track.GetIntLen();
  fBoundaryV[itrack] = track.Boundary();
  fPendingV[itrack] = track.Pending();
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  track.Path()->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  track.NextPath()->CopyTo(fNextpathV[itrack]);
  fNtracks++;
  return itrack;
}

//______________________________________________________________________________
int GeantTrack_v::AddTrackSync(GeantTrack &track) {
  // Add track in a concurrent way. Assumes that this array
  // Is currently being filled while held by the basket manager and NOT being
  // transported.
  // The array has to be compact and should have enough alocated space.
  // Returns the location where the track was added.
  assert(fCompact);
  assert(GetNtracks() < fMaxtracks);
  int itrack = fNtracks++;
  fEventV[itrack] = track.Event();
  fEvslotV[itrack] = track.EventSlot();
  fParticleV[itrack] = track.Particle();
  fMotherV[itrack] = track.Mother();
  fPDGV[itrack] = track.PDG();
  fGVcodeV[itrack] = track.GVcode();
  fEindexV[itrack] = track.EIndex();
  fBindexV[itrack] = track.BIndex();
  fChargeV[itrack] = track.Charge();
  fProcessV[itrack] = track.Process();
  fNstepsV[itrack] = track.GetNsteps();
  fSpeciesV[itrack] = track.Species();
  fStatusV[itrack] = track.Status();
  fMassV[itrack] = track.Mass();
  fXposV[itrack] = track.X();
  fYposV[itrack] = track.Y();
  fZposV[itrack] = track.Z();
  fXdirV[itrack] = track.Dx();
  fYdirV[itrack] = track.Dy();
  fZdirV[itrack] = track.Dz();
  fPV[itrack] = track.P();
  fEV[itrack] = track.E();
  fTimeV[itrack] = track.Time();
  fEdepV[itrack] = track.Edep();
  fPstepV[itrack] = track.GetPstep();
  fStepV[itrack] = track.GetStep();
  fSnextV[itrack] = track.GetSnext();
  fSafetyV[itrack] = track.GetSafety();
  fNintLenV[itrack] = track.GetNintLen();
  fIntLenV[itrack] = track.GetIntLen();
  fBoundaryV[itrack] = track.Boundary();
  fPendingV[itrack] = track.Pending();
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  track.Path()->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  track.NextPath()->CopyTo(fNextpathV[itrack]);
  return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::GetTrack(int i, GeantTrack &track) const {
  // Extract a single track from array.
  track.SetEvent(fEventV[i]);
  track.SetEvslot(fEvslotV[i]);
  track.SetParticle(fParticleV[i]);
  track.SetMother(fMotherV[i]);
  track.SetPDG(fPDGV[i]);
  track.SetGVcode(fGVcodeV[i]);
  track.SetEindex(fEindexV[i]);
  track.SetBindex(fBindexV[i]);
  track.SetCharge(fChargeV[i]);
  track.SetProcess(fProcessV[i]);
  track.SetNsteps(fNstepsV[i]);
  track.SetSpecies(fSpeciesV[i]);
  track.SetStatus(fStatusV[i]);
  track.SetMass(fMassV[i]);
  track.SetPosition(fXposV[i], fYposV[i], fZposV[i]);
  track.SetDirection(fXdirV[i], fYdirV[i], fZdirV[i]);
  track.SetP(fPV[i]);
  track.SetE(fEV[i]);
  track.SetTime(fTimeV[i]);
  track.SetEdep(fEdepV[i]);
  track.SetPstep(fPstepV[i]);
  track.SetStep(fStepV[i]);
  track.SetSnext(fSnextV[i]);
  track.SetSafety(fSafetyV[i]);
  track.SetNintLen(fNintLenV[i]);
  track.SetIntLen(fIntLenV[i]);
  track.SetBoundary(fBoundaryV[i]);
  track.SetPending(fPendingV[i]);
  track.SetPath(fPathV[i]);
  track.SetNextPath(fNextpathV[i]);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::AddTrack(GeantTrack_v &arr, int i, bool /*import*/) {
// Add track from different array
// If addition is done on top of non-compact array,
// the track will be inserted without updating the number of tracks.
// Returns the location where the track was added.
#ifdef VERBOSE
  arr.PrintTrack(i);
#endif
  int itrack = GetNtracks();
  if (!fCompact)
    itrack = fHoles->FirstSetBit();
  if (itrack == fMaxtracks) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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
  fMotherV[itrack] = arr.fMotherV[i];
  fPDGV[itrack] = arr.fPDGV[i];
  fGVcodeV[itrack] = arr.fGVcodeV[i];
  fEindexV[itrack] = arr.fEindexV[i];
  fBindexV[itrack] = arr.fBindexV[i];
  fChargeV[itrack] = arr.fChargeV[i];
  fProcessV[itrack] = arr.fProcessV[i];
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
  fNintLenV[itrack] = arr.fNintLenV[i];
  fIntLenV[itrack] = arr.fIntLenV[i];
  fBoundaryV[itrack] = arr.fBoundaryV[i];
  fPendingV[itrack] = arr.fPendingV[i];
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  arr.fPathV[i]->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  arr.fNextpathV[i]->CopyTo(fNextpathV[itrack]);
  fNtracks++;
  return itrack;
}

//______________________________________________________________________________
int GeantTrack_v::AddTrackSync(GeantTrack_v &arr, int i) {
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
  int itrack = fNtracks++;

  fEventV[itrack] = arr.fEventV[i];
  fEvslotV[itrack] = arr.fEvslotV[i];
  fParticleV[itrack] = arr.fParticleV[i];
  fMotherV[itrack] = arr.fMotherV[i];
  fPDGV[itrack] = arr.fPDGV[i];
  fGVcodeV[itrack] = arr.fGVcodeV[i];
  fEindexV[itrack] = arr.fEindexV[i];
  fBindexV[itrack] = arr.fBindexV[i];
  fChargeV[itrack] = arr.fChargeV[i];
  fProcessV[itrack] = arr.fProcessV[i];
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
  fNintLenV[itrack] = arr.fNintLenV[i];
  fIntLenV[itrack] = arr.fIntLenV[i];
  fBoundaryV[itrack] = arr.fBoundaryV[i];
  fPendingV[itrack] = arr.fPendingV[i];
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  arr.fPathV[i]->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  arr.fNextpathV[i]->CopyTo(fNextpathV[itrack]);
  return itrack;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::AddTrackSyncAt(int itrack, GeantTrack_v &arr, int i) {
  // Add track from different array in a concurrent way. Assumes that this array
  // Is currently being filled while held by the basket manager and NOT being
  // transported.
  // The array has to be compact and should have enough alocated space.
  // Returns the location where the track was added.
  // This *assumes* that no track has been yet added to this array slot.
  // During the use of this routines (likely on a coprocessor where multiple
  // threads are accessing the same GeantTrack_v, the GeantTrack_v is in
  // a somewhat unsual state, where the array is not compact but fHole
  // is not maintain properly (See non concurrent save code comment out
  // below.

  assert(itrack < fMaxtracks);
#ifdef VERBOSE
  arr.PrintTrack(i);
#endif

  // Technically, after setting fHoles to all on,
  // we really should be doing:
  //   fHoles->ResetBitNumber(itrack);
  //   fSelected->ResetBitNumber(itrack);
  // which use bit operation which are not (yet?)
  // done atomically,
  // and we should do:
  //   atomically: fNtracks = max(fNtracks,itrack)

#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  atomicAdd(&fNtracks, 1);
#else
  ++fNtracks;
#endif

  fEventV[itrack] = arr.fEventV[i];
  fEvslotV[itrack] = arr.fEvslotV[i];
  fParticleV[itrack] = arr.fParticleV[i];
  fMotherV[itrack] = arr.fMotherV[i];
  fPDGV[itrack] = arr.fPDGV[i];
  fGVcodeV[itrack] = arr.fGVcodeV[i];
  fEindexV[itrack] = arr.fEindexV[i];
  fBindexV[itrack] = arr.fBindexV[i];
  fChargeV[itrack] = arr.fChargeV[i];
  fProcessV[itrack] = arr.fProcessV[i];
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
  fNintLenV[itrack] = arr.fNintLenV[i];
  fIntLenV[itrack] = arr.fIntLenV[i];
  fBoundaryV[itrack] = arr.fBoundaryV[i];
  fPendingV[itrack] = arr.fPendingV[i];
  // Copy the volume paths
  size_t size_vpath = VolumePath_t::SizeOfInstance(fMaxDepth);
  fPathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + itrack * size_vpath);
  arr.fPathV[i]->CopyTo(fPathV[itrack]);
  fNextpathV[itrack] = reinterpret_cast<VolumePath_t *>(fVPstart + (fMaxtracks + itrack) * size_vpath);
  arr.fNextpathV[i]->CopyTo(fNextpathV[itrack]);
  return itrack;
}

//______________________________________________________________________________
void GeantTrack_v::AddTracks(GeantTrack_v &arr, int istart, int iend, bool /*import*/) {
  // Add tracks from different array. Single thread at a time.
  int ncpy = iend - istart + 1;
  int ntracks = GetNtracks();
  if (ntracks + ncpy >= fMaxtracks) {
    Resize(Math::Max<double>(2 * fMaxtracks, ntracks + ncpy));
  }
  memcpy(&fEventV[ntracks], &arr.fEventV[istart], ncpy * sizeof(int));
  memcpy(&fEvslotV[ntracks], &arr.fEvslotV[istart], ncpy * sizeof(int));
  memcpy(&fParticleV[ntracks], &arr.fParticleV[istart], ncpy * sizeof(int));
  memcpy(&fMotherV[ntracks], &arr.fMotherV[istart], ncpy * sizeof(int));
  memcpy(&fPDGV[ntracks], &arr.fPDGV[istart], ncpy * sizeof(int));
  memcpy(&fGVcodeV[ntracks], &arr.fGVcodeV[istart], ncpy * sizeof(int));
  memcpy(&fEindexV[ntracks], &arr.fEindexV[istart], ncpy * sizeof(int));
  memcpy(&fBindexV[ntracks], &arr.fBindexV[istart], ncpy * sizeof(int));
  memcpy(&fChargeV[ntracks], &arr.fChargeV[istart], ncpy * sizeof(int));
  memcpy(&fProcessV[ntracks], &arr.fProcessV[istart], ncpy * sizeof(int));
  memcpy(&fNstepsV[ntracks], &arr.fNstepsV[istart], ncpy * sizeof(int));
  memcpy(&fSpeciesV[ntracks], &arr.fSpeciesV[istart], ncpy * sizeof(Species_t));
  memcpy(&fStatusV[ntracks], &arr.fStatusV[istart], ncpy * sizeof(TrackStatus_t));
  memcpy(&fMassV[ntracks], &arr.fMassV[istart], ncpy * sizeof(double));
  memcpy(&fXposV[ntracks], &arr.fXposV[istart], ncpy * sizeof(double));
  memcpy(&fYposV[ntracks], &arr.fYposV[istart], ncpy * sizeof(double));
  memcpy(&fZposV[ntracks], &arr.fZposV[istart], ncpy * sizeof(double));
  memcpy(&fXdirV[ntracks], &arr.fXdirV[istart], ncpy * sizeof(double));
  memcpy(&fYdirV[ntracks], &arr.fYdirV[istart], ncpy * sizeof(double));
  memcpy(&fZdirV[ntracks], &arr.fZdirV[istart], ncpy * sizeof(double));
  memcpy(&fPV[ntracks], &arr.fPV[istart], ncpy * sizeof(double));
  memcpy(&fEV[ntracks], &arr.fEV[istart], ncpy * sizeof(double));
  memcpy(&fTimeV[ntracks], &arr.fTimeV[istart], ncpy * sizeof(double));
  memcpy(&fEdepV[ntracks], &arr.fEdepV[istart], ncpy * sizeof(double));
  memcpy(&fPstepV[ntracks], &arr.fPstepV[istart], ncpy * sizeof(double));
  memcpy(&fStepV[ntracks], &arr.fStepV[istart], ncpy * sizeof(double));
  memcpy(&fSnextV[ntracks], &arr.fSnextV[istart], ncpy * sizeof(double));
  memcpy(&fSafetyV[ntracks], &arr.fSafetyV[istart], ncpy * sizeof(double));
  memcpy(&fNintLenV[ntracks], &arr.fNintLenV[istart], ncpy * sizeof(double));
  memcpy(&fIntLenV[ntracks], &arr.fIntLenV[istart], ncpy * sizeof(double));
  memcpy(&fBoundaryV[ntracks], &arr.fBoundaryV[istart], ncpy * sizeof(bool));
  memcpy(&fPendingV[ntracks], &arr.fPendingV[istart], ncpy * sizeof(bool));

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
void GeantTrack_v::SwapTracks(int i, int j) {
  // Swap two tracks in the container
  double tdbl;
  int tint;
  bool tbool;
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
  tint = fMotherV[i];
  fMotherV[i] = fMotherV[j];
  fMotherV[j] = tint;
  tint = fPDGV[i];
  fPDGV[i] = fPDGV[j];
  fPDGV[j] = tint;
  tint = fGVcodeV[i];
  fGVcodeV[i] = fGVcodeV[j];
  fGVcodeV[j] = tint;
  tint = fEindexV[i];
  fEindexV[i] = fEindexV[j];
  fEindexV[j] = tint;
  tint = fBindexV[i];
  fBindexV[i] = fBindexV[j];
  fBindexV[j] = tint;
  tint = fChargeV[i];
  fChargeV[i] = fChargeV[j];
  fChargeV[j] = tint;
  tint = fProcessV[i];
  fProcessV[i] = fProcessV[j];
  fProcessV[j] = tint;
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
  tdbl = fNintLenV[i];
  fNintLenV[i] = fNintLenV[j];
  fNintLenV[j] = tdbl;
  tdbl = fIntLenV[i];
  fIntLenV[i] = fIntLenV[j];
  fIntLenV[j] = tdbl;
  tbool = fBoundaryV[i];
  fBoundaryV[i] = fBoundaryV[j];
  fBoundaryV[j] = tbool;
  tbool = fPendingV[i];
  fPendingV[i] = fPendingV[j];
  fPendingV[j] = tbool;
  tptr = fPathV[i];
  fPathV[i] = fPathV[j];
  fPathV[j] = tptr;
  tptr = fNextpathV[i];
  fNextpathV[i] = fNextpathV[j];
  fNextpathV[j] = tptr;
  bool sel = fSelected->TestBitNumber(j);
  fSelected->SetBitNumber(j, fSelected->TestBitNumber(i));
  fSelected->SetBitNumber(i, sel);
}

//______________________________________________________________________________
void GeantTrack_v::ReplaceTrack(int i, int j) {
  // Replace content of track i with the one of track j
  // WorkloadManager *wm = WorkloadManager::Instance();
  fEventV[i] = fEventV[j];
  fEvslotV[i] = fEvslotV[j];
  fParticleV[i] = fParticleV[j];
  fMotherV[i] = fMotherV[j];
  fPDGV[i] = fPDGV[j];
  fGVcodeV[i] = fGVcodeV[j];
  fEindexV[i] = fEindexV[j];
  fBindexV[i] = fBindexV[j];
  fChargeV[i] = fChargeV[j];
  fProcessV[i] = fProcessV[j];
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
  fNintLenV[i] = fNintLenV[j];
  fIntLenV[i] = fIntLenV[j];
  fBoundaryV[i] = fBoundaryV[j];
  fPendingV[i] = fPendingV[j];
  //   if (!fPathV[i]) fPathV[i] = wm->NavStates()->Borrow();
  //   if (!fNextpathV[i]) fNextpathV[i] = wm->NavStates()->Borrow();
  fPathV[i] = fPathV[j];         // fPathV[j] = 0;
  fNextpathV[i] = fNextpathV[j]; // fNextpathV[j] = 0;
  fSelected->SetBitNumber(i, fSelected->TestBitNumber(j));
}

//______________________________________________________________________________
void GeantTrack_v::DeleteTrack(int /*itr*/) {
  // Delete branch arrays for this track. The track should not have a copy, this has
  // to be called after a killed track is removed by the scheduler.
  //   WorkloadManager *wm = WorkloadManager::Instance();
  //   wm->NavStates()->release(fPathV[itr]);
  //   fPathV[itr] = 0;
  //   wm->NavStates()->release(fNextpathV[itr]);
  //   fNextpathV[itr] = 0;
}

//______________________________________________________________________________
void GeantTrack_v::RemoveTracks(int from, int to) {
// Remove tracks from the container. The method assumes that the tracks were
// copied to another container beforehand.
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (!fCompact)
    Geant::Error("RemoveTracks","%s", "Not compact");
#endif
  int ntracks = GetNtracks();
  if (to >= ntracks - 1) {
    int nzero = ntracks - from;
    memset(&fPathV[from], 0, nzero * sizeof(VolumePath_t *));
    memset(&fNextpathV[from], 0, nzero * sizeof(VolumePath_t *));
  }
  int ncpy = fNtracks - to - 1;
  memmove(&fEventV[from], &fEventV[to + 1], ncpy * sizeof(int));
  memmove(&fEvslotV[from], &fEvslotV[to + 1], ncpy * sizeof(int));
  memmove(&fParticleV[from], &fParticleV[to + 1], ncpy * sizeof(int));
  memmove(&fMotherV[from], &fMotherV[to + 1], ncpy * sizeof(int));
  memmove(&fPDGV[from], &fPDGV[to + 1], ncpy * sizeof(int));
  memmove(&fGVcodeV[from], &fGVcodeV[to + 1], ncpy * sizeof(int));
  memmove(&fEindexV[from], &fEindexV[to + 1], ncpy * sizeof(int));
  memmove(&fBindexV[from], &fBindexV[to + 1], ncpy * sizeof(int));
  memmove(&fChargeV[from], &fChargeV[to + 1], ncpy * sizeof(int));
  memmove(&fProcessV[from], &fProcessV[to + 1], ncpy * sizeof(int));
  memmove(&fNstepsV[from], &fNstepsV[to + 1], ncpy * sizeof(int));
  memmove(&fSpeciesV[from], &fSpeciesV[to + 1], ncpy * sizeof(Species_t));
  memmove(&fStatusV[from], &fStatusV[to + 1], ncpy * sizeof(TrackStatus_t));
  memmove(&fMassV[from], &fMassV[to + 1], ncpy * sizeof(double));
  memmove(&fXposV[from], &fXposV[to + 1], ncpy * sizeof(double));
  memmove(&fYposV[from], &fYposV[to + 1], ncpy * sizeof(double));
  memmove(&fZposV[from], &fZposV[to + 1], ncpy * sizeof(double));
  memmove(&fXdirV[from], &fXdirV[to + 1], ncpy * sizeof(double));
  memmove(&fYdirV[from], &fYdirV[to + 1], ncpy * sizeof(double));
  memmove(&fZdirV[from], &fZdirV[to + 1], ncpy * sizeof(double));
  memmove(&fPV[from], &fPV[to + 1], ncpy * sizeof(double));
  memmove(&fEV[from], &fEV[to + 1], ncpy * sizeof(double));
  memmove(&fTimeV[from], &fTimeV[to + 1], ncpy * sizeof(double));
  memmove(&fEdepV[from], &fEdepV[to + 1], ncpy * sizeof(double));
  memmove(&fPstepV[from], &fPstepV[to + 1], ncpy * sizeof(double));
  memmove(&fStepV[from], &fStepV[to + 1], ncpy * sizeof(double));
  memmove(&fSnextV[from], &fSnextV[to + 1], ncpy * sizeof(double));
  memmove(&fSafetyV[from], &fSafetyV[to + 1], ncpy * sizeof(double));
  memmove(&fNintLenV[from], &fNintLenV[to + 1], ncpy * sizeof(double));
  memmove(&fIntLenV[from], &fIntLenV[to + 1], ncpy * sizeof(double));
  memmove(&fBoundaryV[from], &fBoundaryV[to + 1], ncpy * sizeof(bool));
  memmove(&fPendingV[from], &fPendingV[to + 1], ncpy * sizeof(bool));
  memmove(&fPathV[from], &fPathV[to + 1], ncpy * sizeof(VolumePath_t *));
  memmove(&fNextpathV[from], &fNextpathV[to + 1], ncpy * sizeof(VolumePath_t *));
  fNtracks -= to - from + 1;
  fSelected->ResetAllBits();
  fNselected = 0;
}

//______________________________________________________________________________
int GeantTrack_v::Compact(GeantTrack_v *moveto) {
  // Compact the holes in the array. Return number of active elements. This will
  // lose the track fields in the holes, so information from the holes has to be
  // copied beforehand
  int ntracks = GetNtracks();
  if (ntracks == 0 || fCompact)
    return 0;
  fCompact = true;
  int firsthole = fHoles->FirstSetBit();
  while (firsthole < ntracks) {
    int lastactive = fHoles->LastNullBit(ntracks - 1);
    if (lastactive < ntracks) {
      // move last holes (if any)
      if (moveto && (ntracks - lastactive - 1 > 0))
        moveto->AddTracks(*this, lastactive + 1, ntracks - 1, true);
      ntracks = lastactive + 1;
      if (firsthole == ntracks) {
        SetNtracks(ntracks);
        return ntracks;
      }
    } else {
      // No active tracks left. First copy the hole track to the output
      if (moveto)
        moveto->AddTracks(*this, firsthole, firsthole + ntracks - 1, true);
      SetNtracks(0);
      return 0;
    }
    // replace content of first hole with the last active track
    if (moveto)
      moveto->AddTrack(*this, firsthole, true);
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
int GeantTrack_v::Reshuffle() {
  // Reshuffle tracks according the selection mask. The selected tracks will be
  // moved at the beginning of the array. Tracks should be compacted before.
  if (GetNtracks() == 0)
    return 0;
  fNselected = GetNtracks();
  int firsthole = fSelected->FirstNullBit();
  while (firsthole < fNselected) {
    int lastsel = fSelected->LastSetBit(fNselected - 1);
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
bool GeantTrack_v::Contains(int evstart, int nevents) const {
  // Check if the array contains tracks from a given event range
  int evend = evstart + nevents;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
    if (fEventV[itr] >= evstart && fEventV[itr] < evend)
      return true;
  }
  return false;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::Clear(const char *) {
  // Clear track content and selections
  fNselected = 0;
  int ntracks = GetNtracks();
  if (ntracks) {
    memset(fPathV, 0, ntracks * sizeof(VolumePath_t *));
    memset(fNextpathV, 0, ntracks * sizeof(VolumePath_t *));
  }
  fHoles->ResetAllBits();
  fSelected->ResetAllBits();
  fCompact = true;
  SetNtracks(0);
}

//______________________________________________________________________________
int GeantTrack_v::PropagateStraight(int ntracks, double *crtstep) {
  // Propagate first ntracks along a straight line (neutral particles, no mag.
  // field or for last tiny step). The array has to be reshuffled after selecting
  // the interesting particles using Select method.
  // The crossing tracks get masked as holes in the array.

  // Find next volume
  int icrossed = 0;
  for (int i = 0; i < ntracks; i++) {
    if (fBoundaryV[i]) {
      //*fPathV[i] = *fNextpathV[i];
      fStatusV[i] = kBoundary;
      icrossed++;
    }
  }
  for (int i = 0; i < ntracks; i++) {
    fPstepV[i] -= crtstep[i];
    fSafetyV[i] = 0;
    // Change path to reflect the physical volume for the current track; The
    // relevant path is fPath[i] if the frombdr flag is not set or fNextpath[i]
    // otherwise
    fXposV[i] += crtstep[i] * fXdirV[i];
    fYposV[i] += crtstep[i] * fYdirV[i];
    fZposV[i] += crtstep[i] * fZdirV[i];
#ifdef USE_VECGEOM_NAVIGATOR
//      CheckLocationPathConsistency(i);
#endif
  }
  return icrossed;
}

//______________________________________________________________________________
void GeantTrack_v::PropagateInVolume(int ntracks, const double *crtstep, GeantTaskData *td) {
  // Propagate the selected tracks with crtstep values. The method is to be called
  // only with  charged tracks in magnetic field. The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)
  for (int i = 0; i < ntracks; i++) {
    PropagateInVolumeSingle(i, crtstep[i], td);
  }
}

//______________________________________________________________________________          
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::GetFieldValue( // GeantTaskData *td,
                                 int i, double BfieldOut[3], double *bmagOut) const
{
   vecgeom::Vector3D<double> Position (fXposV[i], fYposV[i], fZposV[i]);
   vecgeom::Vector3D<double> Bfield;
   FieldLookup::GetFieldValue( Position, Bfield, *bmagOut); // , td);
   BfieldOut[0]= Bfield[0];
   BfieldOut[1]= Bfield[1];
   BfieldOut[2]= Bfield[2];   
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::PropagateInVolumeSingle(int i, double crtstep, GeantTaskData * td) {
  // Propagate the selected track with crtstep value. The method is to be called
  // only with  charged tracks in magnetic field.The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)

  // Double_t c = 0.;

  bool useRungeKutta= false;
  GUFieldPropagator *fieldPropagator = nullptr;
  double BfieldInitial[3], bmag= 0.0;
  GetFieldValue(i, BfieldInitial, &bmag);

  useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;
   
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  constexpr bool gPropagator_fUseRK = false; // Temporary work-around until actual implementation ..
  useRungeKutta= gPropagator_fUseRK;   //  Something like this is needed - TBD
#else
  // useRungeKutta= gPropagator->fUseRungeKutta;
  useRungeKutta= td->fPropagator->fConfig->fUseRungeKutta;  
  if( useRungeKutta ){
    fieldPropagator = td->fFieldPropagator;
    assert( fieldPropagator );
  }
#endif

  // Reset relevant variables - TBC: check if changes are needed after the endpoint is knownn
  fStatusV[i] = kInFlight;
  fPstepV[i] -= crtstep;
  if (fPstepV[i] < 1.E-10) {
    fPstepV[i] = 0;
    fStatusV[i] = kPhysics;
  }
  fSafetyV[i] -= crtstep;
  if (fSafetyV[i] < 1.E-10)
    fSafetyV[i] = 0;
  fSnextV[i] -= crtstep;
  if (fSnextV[i] < 1.E-10) {
    fSnextV[i] = 0;
    if (fBoundaryV[i]) {
      fStatusV[i] = kBoundary;
    }
  }
  fStepV[i] += crtstep;

  double curvaturePlus= fabs(GeantTrack::kB2C * fChargeV[i] * bmag) / (fPV[i] + 1.0e-30);  // norm for step
  // 'Curvature' along the full track - not just in the plane perpendicular to the B-field vector

  constexpr double numRadiansMax= 10.0;   //  Too large an angle - many RK steps.  Potential change -> 2.0*PI;
  constexpr double numRadiansMin= 0.05;   //  Very small an angle - helix is adequate.  TBC: Use average B-field value?
      //  A track turning more than 10 radians will be treated approximately
  const double angle= crtstep * curvaturePlus;
  bool mediumAngle = ( numRadiansMin < angle ) && ( angle < numRadiansMax );
  useRungeKutta = useRungeKutta && (mediumAngle);

  bool dominantBz =  std::fabs( std::fabs(BfieldInitial[2]) )
         > 1.e3 *
     std::max( std::fabs( BfieldInitial[0]), std::fabs(BfieldInitial[1]) );

#ifdef DEBUG_FIELD
  printf("--PropagateInVolumeSingle: \n");
  printf("Curvature= %8.4g   CurvPlus= %8.4g  step= %f   Bmag=%8.4g   momentum mag=%f  angle= %g\n"
         Curvature(td, i), curvaturePlus, crtstep, bmag, fPV[i], angle );
#endif
  
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
  using ThreeVector = vecgeom::Vector3D<double>;
  // typedef vecgeom::Vector3D<double>  ThreeVector;
  ThreeVector Position(fXposV[i], fYposV[i], fZposV[i]);
  ThreeVector Direction(fXdirV[i], fYdirV[i], fZdirV[i]);
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

  int propagationType= 0;

  if( useRungeKutta )
  {
     // crtstep = 1.0e-4;   printf( "Setting crtstep = %f -- for DEBUGing ONLY.", crtstep );
     propagationType= 1;
     // PrintTrack(i);
#ifndef VECCORE_CUDA
     fieldPropagator->DoStep(Position,    Direction,    fChargeV[i], fPV[i], crtstep,
                             PositionNew, DirectionNew); // BfieldInitial );

     const bool fCheckingStep= false;
     if( fCheckingStep ) {
        const double epsDiff= 2.0e-3; // bool verbDiff= true );
        StepChecker EndChecker( epsDiff, epsDiff * crtstep, true );
        vecgeom::Vector3D<double> Bfield( BfieldInitial[0], BfieldInitial[1], BfieldInitial[2] );
        EndChecker.CheckStep( Position, Direction, fChargeV[i], fPV[i], crtstep,
                              PositionNew, DirectionNew, Bfield );
     }
     // CheckDirection(DirectionNew);
#endif
  } else {
     constexpr double toKiloGauss= 1.0e+14; // Converts to kilogauss -- i.e. 1 / Unit::kilogauss
                                            // Must agree with values in magneticfield/inc/Units.h
     double Bz = BfieldInitial[2] * toKiloGauss;
     if ( dominantBz ) { // Oldest - constant field along z        
        // printf("h"); std::cout << "h";
        propagationType= 2;
        // Printf("Called Helix-Bz.  Bz= %g , ( Bx = %g, By= %g ) Kilo Gauss\n", Bz, Bx, By );

        // Constant field in Z-direction
        ConstBzFieldHelixStepper stepper( Bz ); //
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  fChargeV[i], fPV[i], crtstep,
                                               PositionNew, DirectionNew);
     } else {
        propagationType= 3;
        // double Bx = BfieldInitial[0] * toKiloGauss;
        // double By = BfieldInitial[1] * toKiloGauss;
        // printf("H"); std::cout << "H";
        if( ! CheckDirection( i, 1.0e-4 ) )
           PrintTrack(i, "Failed check of *direction* - input to General Helix stepper.");

        // Printf("Called Helix-General.  Bz= %g , Bx = %g, By= %g ", Bz, Bx, By );
        
        Geant::ConstFieldHelixStepper stepper( BfieldInitial );
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  fChargeV[i], fPV[i], crtstep,
                                               PositionNew, DirectionNew);
     }
  }

  fXposV[i] = PositionNew.x();
  fYposV[i] = PositionNew.y();
  fZposV[i] = PositionNew.z();

  // double oldMag = DirectionNew.Mag();
  // if( std::fabs( oldMag - 1.0 ) > 1.e-6 ) { fNumBadNormals++; fMaxBadNormalFld= std::max(fMaxBadNormal, oldMag); fMinBadNormalFld = std::min(fMinBadNormal, oldMag); }
  
  // Normalize direction here - to avoid surprises  ( but it limits ability to check bad integration )
  ThreeVector DirectionUnit = DirectionNew.Unit();
  fXdirV[i] = DirectionUnit.x();
  fYdirV[i] = DirectionUnit.y();
  fZdirV[i] = DirectionUnit.z();

#ifdef REPORT_AND_CHECK
  double newMag= DirectionNew.Mag();
  Printf(" -- State after propagation in field:  Position= %f, %f, %f   Direction= %f, %f, %f  - mag old, new, new-1.0 = %10.8f %10.8f %7.2g",
         fXposV[i], fYposV[i], fZposV[i],
         fXdirV[i], fYdirV[i], fZdirV[i], oldMag, newMag, newMag-1.0 );
         // DirectionNew.Mag()-1.0  );

  const char* Msg[4]= { "After propagation in field - type Unknown(ERROR) ",
                        "After propagation in field - with RK           ",
                        "After propagation in field - with Helix-Bz     ",
                        "After propagation in field - with Helix-General" };
  CheckTrack(i, Msg[propagationType] );
#endif

#if 0
  ThreeVector SimplePosition = Position + crtstep * Direction;
  // double diffpos2 = (PositionNew - Position).Mag2();
  double diffpos2 = (PositionNew - SimplePosition).Mag2();
  //   -- if (Math::Sqrt(diffpos)>0.01*crtstep) {
  const double drift= 0.01*crtstep;
  if ( diffpos2>drift*drift ){
r      double diffpos= Math::Sqrt(diffpos2);
      // Geant::Print("PropagateInVolumeSingle","relative difference in pos = %g", diffpos/crtstep);
      Geant::Print("PropagateInVolumeSingle","difference in pos = %g (abs) %g (relative) , step= %g",
                   diffpos, diffpos/crtstep, crtstep);
  }
#endif
}

#ifdef USE_VECGEOM_NAVIGATOR
void GeantTrack_v::CheckLocationPathConsistency(int itr) const {
  VECGEOM_NAMESPACE::NavigationState *a =
      VECGEOM_NAMESPACE::NavigationState::MakeInstance(VECGEOM_NAMESPACE::GeoManager::Instance().getMaxDepth());
  a->Clear();
  VECGEOM_NAMESPACE::SimpleNavigator nav;
  nav.LocatePoint(VECGEOM_NAMESPACE::GeoManager::Instance().GetWorld(),
                  VECGEOM_NAMESPACE::Vector3D<VECGEOM_NAMESPACE::Precision>(fXposV[itr], fYposV[itr], fZposV[itr]), *a,
                  true);
  if (a->Top() != NULL && a->Top() != fPathV[itr]->Top()) {
    Geant::Print("","INCONSISTENT LOCATION PATH PAIR PRODUCED FOR TRACK %d", itr);
#ifdef VECGEOM_ROOT
    Geant::Print("","REAL");
    a->GetCurrentNode()->Print();
    Geant::Print("","REPORTED");
    fPathV[itr]->GetCurrentNode()->Print();
//  printrace();
#endif
  }

  // release object
  VECGEOM_NAMESPACE::NavigationState::ReleaseInstance(a);
}
#endif

//______________________________________________________________________________
int GeantTrack_v::SortByStatus(TrackStatus_t status) {
  // Sort tracks by a given status.
  int nsel = 0;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
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
int GeantTrack_v::SortByLimitingDiscreteProcess() {
  // Sort tracks for which the step was limited by discrete processes.
  int nsel = 0;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
    if (fStatusV[itr] == kPhysics && fEindexV[itr] == 1000) {
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
int GeantTrack_v::RemoveByStatus(TrackStatus_t status, GeantTrack_v &output) {
  // Remove tracks with given status from the container to the output vector,
  // then compact.
  int nremoved = 0;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
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
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::PrintTrack(int itr, const char *msg) const {
  // Print info for a given track
  const char *status[8] = {"alive", "killed", "inflight", "boundary", "exitSetup", "physics", "postponed", "new"};
#ifdef USE_VECGEOM_NAVIGATOR
  Geant::Print(msg,
      "== Track %d: evt=%d slt=%d part=%d mth=%d pdg=%d gVc=%d chg=%d proc=%d nstp=%d spc=%d status=%s mass=%g "
      "xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n",
      itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fMotherV[itr], fPDGV[itr], fGVcodeV[itr],
      fChargeV[itr], fProcessV[itr], fNstepsV[itr], (int)fSpeciesV[itr], status[int(fStatusV[itr])],
      fMassV[itr], fXposV[itr], fYposV[itr], fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr],
      fTimeV[itr], fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fNintLenV[itr], fIntLenV[itr], fBoundaryV[itr]);

#ifndef VECCORE_CUDA
  fPathV[itr]->Print();
  fNextpathV[itr]->Print();
#endif
#else
  TString path;
  fPathV[itr]->GetPath(path);
  TString nextpath;
  fNextpathV[itr]->GetPath(nextpath);

  Geant::Print(msg, "== Track %d: evt=%d slt=%d part=%d mth=%d pdg=%d gVc=%d eind=%d chg=%d proc=%d nstp=%d "
         "spc=%d status=%s mass=%g xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g "
         "time=%g edep=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n pth=%s npth=%s\n",
         itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fMotherV[itr], fPDGV[itr], fEindexV[itr], fGVcodeV[itr],
         fChargeV[itr], fProcessV[itr], fNstepsV[itr], (int)fSpeciesV[itr], status[int(fStatusV[itr])],
         fMassV[itr], fXposV[itr], fYposV[itr], fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr],
         fTimeV[itr], fEdepV[itr], fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fNintLenV[itr], fIntLenV[itr], fBoundaryV[itr], path.Data(),
         nextpath.Data());
#endif
}

//______________________________________________________________________________
void GeantTrack_v::PrintTracks(const char *msg) const {
  // Print all tracks
  int ntracks = GetNtracks();
  Geant::Print(msg,"%s","");
  for (int i = 0; i < ntracks; i++)
    PrintTrack(i);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::ComputeTransportLength(int ntracks, GeantTaskData *td) {
// Vector version for proposing the geometry step. All tracks have to be in
// the same volume
#ifdef USE_VECGEOM_NAVIGATOR
//#define VECTORIZED_GEOMETRY
#ifdef VECTORIZED_GEOMETRY
  // We reshuffle tracks for which the current safety allows for the proposed step
  int nsel = 0;
  for (int itr = 0; itr < ntracks; ++itr) {
    // Check if current safety allows for the proposed step
    if (fSafetyV[itr] > fPstepV[itr]) {
      fSnextV[itr] = fPstepV[itr];
      fBoundaryV[itr] = false;
      *fNextpathV[itr] = *fPathV[itr];
      continue;
    } else {
      // These tracks have to be processed vectorized by geometry
      Select(itr);
      nsel++;
    }
  }
  if (nsel == ntracks) {
    // All tracks have proposed steps bigger than safety -> unmark them
    DeselectAll();
  } else {
    if (nsel>0) Reshuffle();
    else return;
  }
  // Printf("In vec find next boundary and step\n");

  // call the vector interface
  VectorNavInterface
    ::NavFindNextBoundaryAndStep(nsel, fPstepV,
                              fXposV, fYposV, fZposV,
                              fXdirV, fYdirV, fZdirV,
                              (const VolumePath_t **)fPathV, fNextpathV,
                              fSnextV, fSafetyV, fBoundaryV);

   // Update number of calls to geometry
   td->fNsnext += 1;
#else
  // Non-vectorized looped implementation
#ifdef NEW_NAVIGATION
  ScalarNavInterfaceVGM
#else
  ScalarNavInterfaceVG
#endif // NEW_NAVIGATION
   ::NavFindNextBoundaryAndStep(ntracks, fPstepV,
                              fXposV,fYposV,fZposV,
                              fXdirV, fYdirV, fZdirV,
                              (const VolumePath_t**)fPathV, fNextpathV,
                              fSnextV, fSafetyV, fBoundaryV);
  // Update number of calls to geometry (consider N scalar calls)
  td->fNsnext += ntracks;
#endif // VECTORIZED_GEOMETRY
  // perform a couple of additional checks/ set status flags and so on
  for (int itr = 0; itr < ntracks; ++itr) {
    if ((fNextpathV[itr]->IsOutside() && fSnextV[itr] < 1.E-6) || fSnextV[itr] > 1.E19)
      fStatusV[itr] = kExitingSetup;
  }
#else
  // TGeo implementation fall on looped version
  for (int i = 0; i < ntracks; ++i) {
    ComputeTransportLengthSingle(i, td);
  }
#endif // USE_VECGEOM_NAVIGATOR
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::ComputeTransportLengthSingle(int itr, GeantTaskData *td) {
// Computes snext and safety for a single track. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.

#ifdef USE_VECGEOM_NAVIGATOR
//#define NEW_NAVIGATION
#ifdef NEW_NAVIGATION
  ScalarNavInterfaceVGM
   ::NavFindNextBoundaryAndStep(1, &fPstepV[itr],&fXposV[itr],&fYposV[itr],&fZposV[itr],
                                &fXdirV[itr],&fYdirV[itr],&fZdirV[itr],
                                (const VolumePath_t**)(&fPathV[itr]), &fNextpathV[itr],
                                &fSnextV[itr],&fSafetyV[itr],&fBoundaryV[itr]);
#else
  ScalarNavInterfaceVG
   ::NavFindNextBoundaryAndStep(1, &fPstepV[itr],&fXposV[itr],&fYposV[itr],&fZposV[itr],
                                &fXdirV[itr],&fYdirV[itr],&fZdirV[itr],
                                (const VolumePath_t**)(&fPathV[itr]), &fNextpathV[itr],
                                &fSnextV[itr],&fSafetyV[itr],&fBoundaryV[itr]);
#endif // NEW_NAVIGATION
#else
// ROOT geometry
  ScalarNavInterfaceTGeo
   ::NavFindNextBoundaryAndStep(1, &fPstepV[itr],&fXposV[itr],&fYposV[itr],&fZposV[itr],
                                &fXdirV[itr],&fYdirV[itr],&fZdirV[itr],
                                (const VolumePath_t**)(&fPathV[itr]), &fNextpathV[itr],
                                &fSnextV[itr],&fSafetyV[itr],&fBoundaryV[itr]);
#endif // USE_VECGEOM_NAVIGATOR
  // Update number of calls to geometry
  td->fNsnext++;
  // if outside detector or enormous step mark particle as exiting the detector
  if (fNextpathV[itr]->IsOutside() || fSnextV[itr] > 1.E19)
    fStatusV[itr] = kExitingSetup;
}

//______________________________________________________________________________
TransportAction_t GeantTrack_v::PostponedAction(int ntracks) const {
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
int GeantTrack_v::PropagateTracks(GeantTaskData *td) {
  // Propagate the ntracks in the current volume with their physics steps (already
  // computed)
  // Vectors are pushed downstream when efficient.
  GeantTrack_v &output = *td->fTransported;
  int ntracks = GetNtracks();
  // Check if tracking the remaining tracks can be postponed
  TransportAction_t action = PostponedAction(ntracks);
  if (action == kPostpone) {
    PostponeTracks(output);
    return 0;
  }
  if (action != kVector)
    return PropagateTracksScalar(td, 0);
// Compute transport length in geometry, limited by the physics step
#ifdef BUG_HUNT
  GeantPropagator *prop = td->fPropagator;
  BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "PropagateTracks");
#endif
  ComputeTransportLength(ntracks, td);
  //     Printf("====== After ComputeTransportLength:");
  //     PrintTracks();
  // double sumEin=0.0, sumEout=0.0; // , sumEdep=0.0
  // for (int ix = 0; ix < ntracks; ix++) { sumEin += fEV[ix]; }

#ifdef BUG_HUNT
  BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterCompTransLen");
#endif

  int itr = 0;
  int icrossed = 0;
  int nsel = 0;
  double lmax;
  const double eps = 1.E-2; // 100 micron
  // const double bmag = prop->fConfig->fBmag;
  double Bfield[3], bmag= 0.0;

  unsigned int numNeutral= 0, numCharged=0, numStraight=0, numCurved=0; // , numPhysics=0;

  // Remove dead tracks, propagate neutrals
  for (itr = 0; itr < ntracks; itr++) {
    // Mark dead tracks for copy/removal
    if (fSnextV[itr] < 0) {
      Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", fParticleV[itr]);
      PrintTrack(itr);
      fStatusV[itr] = kKilled;
    }
    if (fStatusV[itr] == kKilled) {
      MarkRemoved(itr);
      continue;
    }
    // Propagate straight tracks to the precomputed location and update state,
    // then mark them for copy/removal
    // (Inlined from PropagateStraight)

    // if (fChargeV[itr] == 0 || bmag < 1.E-10) {
    
    bool straightTraj= ( fChargeV[itr] == 0 );
    if( !straightTraj ) {
       numCharged++;
       GetFieldValue( itr, Bfield, &bmag);
       // td->StoreFieldValue(itr, Bfield, bmag);   // Store it in Task-Data array !?
       straightTraj = bmag < 1.E-10 * fieldUnits::kilogauss;
       // printf("bmag = %9.3g kiloGauss\n", bmag / fieldUnits::kilogauss );
    } else {
       // td->ClearFieldValue(itr);
       numNeutral++;
    }
    if( straightTraj ) {
      numStraight++;
       
      // Do straight propagation to physics process or boundary
      if (fBoundaryV[itr]) {
        if (fNextpathV[itr]->IsOutside())
          fStatusV[itr] = kExitingSetup;
        else
          fStatusV[itr] = kBoundary;
        icrossed++;
      } else {
        fStatusV[itr] = kPhysics;
        // Update number of steps to physics
        td->fNphys++;
      }
      fPstepV[itr] -= fSnextV[itr];
      fStepV[itr] += fSnextV[itr];
      fSafetyV[itr] -= fSnextV[itr];
      if (fSafetyV[itr] < 0.)
        fSafetyV[itr] = 0;
      fXposV[itr] += fSnextV[itr] * fXdirV[itr];
      fYposV[itr] += fSnextV[itr] * fYdirV[itr];
      fZposV[itr] += fSnextV[itr] * fZdirV[itr];
      // Update total number of steps
      td->fNsteps++;
      if (fSnextV[itr] < 1.E-8) td->fNsmall++;
      MarkRemoved(itr);
      fSnextV[itr] = 0;
#ifdef USE_VECGEOM_NAVIGATOR
//            CheckLocationPathConsistency(itr);
#endif
    } else {
       numCurved++;
    }
  }

  /**
  for (int ix = 0; ix < ntracks; ix++) { sumEout += fEV[ix]; }
  for (int ix = 0; ix < output.GetNtracks(); ix++) { sumEout += output.fEV[ix]; }
  if( sumEout - sumEin > 0.001 * sumEin )
     Printf("PropagateTracks: Ein= %8.3g                      Eout= %8.3g                      Balance= %8.3g",
            sumEin, sumEout, sumEout - sumEin );
   **/
  
  // Compact remaining tracks and move the removed oned to the output container
  if (!fCompact)
    Compact(&output);

  // Check if tracking the remaining tracks can be postponed
  action = PostponedAction(fNtracks);

  // static unsigned long totalTracks=0, totalNeutral=0, totalCharged=0, totalStraight=0, totalCurved=0, totalPhysics=0;

  /*****
  unsigned int numCalls=0;
  const unsigned modCalls = 1000;
  if( (++numCalls) % modCalls == 0 ) {
     Printf("\nPropagateTracks: # tracks: Neutral=%4d, Charged=%4d, numStraight=%4d, numCurved=%4d, numPhysics=%4d . Action= %2d",
            numNeutral, numCharged, numStraight, numCurved, numPhysics, action );
  }
  totalTracks += numCharged + numNeutral;
  totalNeutral += numNeutral;
  totalCharged += numCharged;
  totalStraight += numStraight;
  totalCurved  += numCurved;
  totalPhysics += numPhysics;
   *****/
  
  switch (action) {
  case kDone:
    return icrossed;
  case kSingle:
    icrossed += PropagateTracksScalar(td, 1);
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
  double *steps = td->GetDblArray(ntracks);
  for (itr = 0; itr < fNtracks; itr++) {
    lmax = SafeLength(itr, eps);
    lmax = Math::Max<double>(lmax, fSafetyV[itr]);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] =
        (fBoundaryV[itr]) ? Math::Min<double>(lmax, Math::Max<double>(fSnextV[itr], 1.E-4)) : Math::Min<double>(lmax, fPstepV[itr]);
    //    if (fBoundaryV[itr] && steps[itr]<1.E-8) steps[itr] = 1.E-3;
    // Printf("snext=%g lmax=%g", fSnextV[itr], lmax);
    //      Printf("track %d: step=%g (safelen=%g)", itr, steps[itr], lmax);
  }
  // Propagate the vector of tracks
  PropagateInVolume(ntracks, steps, td);
  //Update number of partial steps propagated in field
  td->fNmag += ntracks;
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
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
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
    // Check if boundary have ben crossed
    bool *same = td->GetBoolArray(nsel);
#ifdef USE_VECGEOM_NAVIGATOR
    NavigationState *tmpstate = td->GetPath();
#ifdef VECTORIZED_GEOMETRY
    VectorNavInterface
#else
#ifdef NEW_NAVIGATION
    ScalarNavInterfaceVGM
#else
// Old navigation system
    ScalarNavInterfaceVG
#endif // NEW_NAVIGATION
#endif // VECTORIZED_GEOMETRY
      ::NavIsSameLocation(ntracks, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV, (const VolumePath_t**)fPathV, fNextpathV, same, tmpstate);
#else
// ROOT navigation
    ScalarNavInterfaceTGeo
      ::NavIsSameLocation(ntracks, fXposV, fYposV, fZposV, fXdirV, fYdirV, fZdirV, (const VolumePath_t**)fPathV, fNextpathV, same);
#endif // USE_VECGEOM_NAVIGATOR
    for (itr = 0; itr < nsel; itr++) {
      if (same[itr]) {
        fBoundaryV[itr] = false;
        continue;
      }
      // Boundary crossed
      fStatusV[itr] = kBoundary;
      if (fNextpathV[itr]->IsOutside())
        fStatusV[itr] = kExitingSetup;
      fBoundaryV[itr] = true;
      icrossed++;
      // Update number of steps to physics and total number of steps
      td->fNsteps++;
      if (fStepV[itr] < 1.E-8) td->fNsmall++;
      MarkRemoved(itr);
    }
    //         Printf("====== After finding crossing tracks (ncross=%d):", icrossed);
    //         PrintTracks();
    if (!fCompact)
      Compact(&output);
  }
#ifdef BUG_HUNT
  BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateTracks");
#endif
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::PropagateSingleTrack(int itr, GeantTaskData *td, int stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  // printf(" PropagateSingleTrack called for itr= %d -- Track: ", itr );
  // PrintTrack(itr);
  // CheckTrack(itr, " PropagateSingleTrack called.");
  // PrintTrack(itr);

  int icrossed = 0;
  double step, lmax;
  const double eps = 1.E-2; // 1 micron
  // Compute transport length in geometry, limited by the physics step
  double Bfield[3], bmag= 0.0;
#ifdef BUG_HUNT
  GeantPropagator *prop = td->fPropagator;
  BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "PropagateSingle", itr);
#endif
  ComputeTransportLengthSingle(itr, td);
#ifdef BUG_HUNT
  BreakOnStep(0, 15352, 0, 10, "AfterCompTranspLenSingle");
#endif
  // Mark dead tracks for copy/removal
  if (fSnextV[itr] < 0) {
    Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", fParticleV[itr]);
    PrintTrack(itr);
    fStatusV[itr] = kKilled;
  }
  if (fStatusV[itr] == kKilled) {
    MarkRemoved(itr);
    return icrossed;
  }
  // Stage 0: straight propagation
  if (stage == 0) {
    bool neutral = (fChargeV[itr] == 0);
    if( !neutral ) {
       // printf( " PropagateSingleTrack> getting Field. Charge= %3d ", fChargeV[itr]);
       GetFieldValue(itr, Bfield, &bmag);
       // if( bmag < 1.E-10) { printf(" Tiny field - mag = %g at %f %f %f\n",
       //                             bmag,  fXposV[itr],  fYposV[itr],  fZposV[itr]); }
    }
    // if (fChargeV[itr] == 0 || bmag < 1.E-10) {
    if ( neutral ) { // || bmag < 1.E-10 * kiloGauss ) {
      // Do straight propagation to physics process or boundary
      if (fBoundaryV[itr]) {
        //*fPathV[itr] = *fNextpathV[itr];
        if (fNextpathV[itr]->IsOutside())
          fStatusV[itr] = kExitingSetup;
        else
          fStatusV[itr] = kBoundary;
        icrossed++;
      } else {
        fStatusV[itr] = kPhysics;
        // Update number of steps to physics
        td->fNphys++;
      }
      fPstepV[itr] -= fSnextV[itr];
      fStepV[itr] += fSnextV[itr];
      fSafetyV[itr] -= fSnextV[itr];
      if (fSafetyV[itr] < 0.)
        fSafetyV[itr] = 0;
      fXposV[itr] += fSnextV[itr] * fXdirV[itr];
      fYposV[itr] += fSnextV[itr] * fYdirV[itr];
      fZposV[itr] += fSnextV[itr] * fZdirV[itr];
      // Update total number of steps
      td->fNsteps++;
      if (fSnextV[itr] < 1.E-8) td->fNsmall++;
      fSnextV[itr] = 0;
      MarkRemoved(itr);
#ifdef USE_VECGEOM_NAVIGATOR
//            CheckLocationPathConsistency(itr);
#endif
#ifdef BUG_HUNT
      BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateSingleNeutral",
                  itr);
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
    lmax = Math::Max<double>(lmax, fSafetyV[itr]);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if frombdr=1) or the proposed  physics
    // step (frombdr=0)
    step = (fBoundaryV[itr]) ? Math::Min<double>(lmax, Math::Max<double>(fSnextV[itr], 1.E-4)) :
                               Math::Min<double>(lmax, fPstepV[itr]);
    //      Printf("track %d: step=%g (safelen=%g)", itr, step, lmax);
    // int stepNum= fNstepsV[itr];
    // Printf("track %d: Step #=%3d len=%g proposed=%g (safelen=%9.3g) bndrFlg= %d distLin=%g  ",
    //    itr, stepNum, step, fPstepV[itr], lmax, fFrombdrV[itr], fSnextV[itr] );
    //  Printf("track %d: step=%g (safelen=%g)", itr, step, lmax);
    PropagateInVolumeSingle(itr, step, td);
    //Update number of partial steps propagated in field
    td->fNmag++;
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
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      return icrossed;
    }
    // Select tracks that are in flight or were propagated to boundary with
    // steps bigger than safety
    if (fSafetyV[itr] < 1.E-10 || fSnextV[itr] < 1.E-10) {
      // Check if boundary has been crossed
      bool same = true;
#ifdef USE_VECGEOM_NAVIGATOR
      NavigationState *tmpstate = td->GetPath();
#ifdef VECTORIZED_GEOMETRY
      VectorNavInterface
#else
#ifdef NEW_NAVIGATION
      ScalarNavInterfaceVGM
#else
// Old navigation system
      ScalarNavInterfaceVG
#endif // NEW_NAVIGATION
#endif // VECTORIZED_GEOMETRY
        ::NavIsSameLocation(1, &fXposV[itr], &fYposV[itr], &fZposV[itr], &fXdirV[itr], &fYdirV[itr], &fZdirV[itr], (const VolumePath_t**)(&fPathV[itr]), &fNextpathV[itr], &same, tmpstate);
#else
// ROOT navigation
      ScalarNavInterfaceTGeo
        ::NavIsSameLocation(1, &fXposV[itr], &fYposV[itr], &fZposV[itr], &fXdirV[itr], &fYdirV[itr], &fZdirV[itr], (const VolumePath_t**)(&fPathV[itr]), &fNextpathV[itr], &same);
#endif // USE_VECGEOM_NAVIGATOR
      if (same) {
        fBoundaryV[itr] = false;
        return icrossed;
      }
      // Boundary crossed
      fStatusV[itr] = kBoundary;
      if (fNextpathV[itr]->IsOutside())
        fStatusV[itr] = kExitingSetup;
      fBoundaryV[itr] = true;
      icrossed++;
      // Update number of steps to physics and total number of steps
      td->fNsteps++;
      if (fStepV[itr] < 1.E-8) td->fNsmall++;
      MarkRemoved(itr);
    }
#ifdef USE_VECGEOM_NAVIGATOR
//         CheckLocationPathConsistency(itr);
#endif
  }
#ifdef BUG_HUNT
  BreakOnStep(prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateSingle", itr);
#endif
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::PropagateTracksScalar(GeantTaskData *td, int stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  GeantTrack_v &output = *td->fTransported;
#endif
  int icrossed = 0;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
    icrossed += PropagateSingleTrack(itr, td, stage);
  }
//   Printf("====== After finding crossing tracks (ncross=%d):", icrossed);
//   PrintTracks();
// Compact remaining tracks and move the removed oned to the output container
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  if (!fCompact)
    Compact(&output);
#endif
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
double GeantTrack_v::Curvature(int i) const {
  using ThreeVector_d = vecgeom::Vector3D<double>;
  
  // Curvature for general field
  const double tiny = 1.E-30;

  double Bfield[3], bmag= 0.0;
  
  GetFieldValue(i, Bfield, &bmag);
  ThreeVector_d MagFld(Bfield[0], Bfield[1], Bfield[2]);

  //  Calculate transverse momentum 'Pt' for field 'B'
  // 
  ThreeVector_d Momentum( fXdirV[i], fXdirV[i], fXdirV[i] );
  Momentum *= fPV[i];
  ThreeVector_d PtransB;  //  Transverse wrt direction of B
  double ratio = 0.0;
  if( bmag > 0 ) ratio = Momentum.Dot( MagFld ) / (bmag*bmag);
  PtransB = Momentum - ratio * MagFld ;
  double Pt_mag = PtransB.Mag();

  return fabs(GeantTrack::kB2C * fChargeV[i] * bmag / (Pt_mag + tiny));
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
double GeantTrack_v::SafeLength(int i, double eps) {
  // Returns the propagation length in field such that the propagated point is
  // shifted less than eps with respect to the linear propagation.
  double c = Curvature(i);
  if (c < 1.E-10)
    return 1.E50;
  return 2. * sqrt(eps / c);
}

//______________________________________________________________________________
int GeantTrack_v::PostponeTracks(GeantTrack_v &output) {
  // Postpone transport of remaining tracks and copy them to the output.
  int npostponed = GetNtracks();
  for (int itr = 0; itr < npostponed; itr++)
    fStatusV[itr] = kPostponed;
  // Move these tracks to the output container
  output.AddTracks(*this, 0, npostponed - 1, true);
  RemoveTracks(0, npostponed - 1);
  Clear();
  return npostponed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeantTrack_v::PostponeTrack(int itr, GeantTrack_v &output) {
  // Postpone transport of a track and copy it to the output.
  // Returns where in the output the track was added.

  fStatusV[itr] = kPostponed;
  // Move these tracks to the output container
  int new_itr = output.AddTrack(*this, itr, true);
  MarkRemoved(itr);
  return new_itr;
}


//______________________________________________________________________________
Volume_t const*GeantTrack_v::GetNextVolume(int i) const {
  // Next volume the track is getting into
#ifdef USE_VECGEOM_NAVIGATOR
  return fNextpathV[i]->Top()->GetLogicalVolume();
#else
  return fNextpathV[i]->GetCurrentNode()->GetVolume();
#endif
}

//______________________________________________________________________________
Volume_t const*GeantTrack_v::GetVolume(int i) const {
  // Current volume the track is into
#ifdef USE_VECGEOM_NAVIGATOR
  return fPathV[i]->Top()->GetLogicalVolume();
#else
  return (fPathV[i]->GetCurrentNode()->GetVolume());
#endif
}

//______________________________________________________________________________
Material_t *GeantTrack_v::GetMaterial(int i) const {
  // Current material the track is into
#ifdef USE_VECGEOM_NAVIGATOR
  Material_t *mat = (Material_t *)GetVolume(i)->GetMaterialPtr();
  return mat;
#else
  Medium_t *med = (Medium_t *)GetVolume(i)->GetMedium();
  // TODO: better to use assert
  if (!med)
    return nullptr;
  return med->GetMaterial();
#endif
}

//______________________________________________________________________________
#ifdef USE_VECGEOM_NAVIGATOR
bool GeantTrack_v::CheckNavConsistency(int /*itr*/) {
  // TO IMPLEMENT WITH VECGEOM
#else
bool GeantTrack_v::CheckNavConsistency(int itr) {
// Check consistency of navigation state for a given track.
// Debugging purpose
  double point[3], local[3];
  point[0] = fXposV[itr];
  point[1] = fYposV[itr];
  point[2] = fZposV[itr];
  fPathV[itr]->GetMatrix()->MasterToLocal(point, local);
  TGeoShape *shape = fPathV[itr]->GetCurrentNode()->GetVolume()->GetShape();
  int evt = fEventV[itr];
  int trk = fParticleV[itr];
  int stp = fNstepsV[itr];
  bool onbound = fBoundaryV[itr];
  bool inside = shape->Contains(local);
  double safin = shape->Safety(local, true);
  double safout = shape->Safety(local, false);

  // 1. Check that the current state really contains the particle position if the position is not declared on boundary.
  if (!onbound && !inside && safout > 0.01) {
    Printf("ERRINSIDE: evt=%d trk=%d stp=%d (%16.14f,  %16.14f, %16.14f) not inside. safout=%g", evt, trk, stp,
           point[0], point[1], point[2], safout);
    //    PrintTrack(itr);
    return false;
  }
  // 2. Check that the safety state is consistent
  if (!onbound && inside) {
    if (safin < fSafetyV[itr] - 1.E-8) {
      Printf("ERRSAFIN: evt=%d trk=%d stp=%d (%16.14f,  %16.14f, %16.14f) safin=%g smaller than track safety=%g", evt,
             trk, stp, point[0], point[1], point[2], safin, fSafetyV[itr]);
      return false;
    }
  }
#endif
  return true;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool GeantTrack_v::BreakOnStep(int evt, int trk, int stp, int nsteps, const char *msg, int itr) {
  // Return true if container has a track with a given number doing a given step from a given event
  // Debugging purpose
  int ntracks = GetNtracks();
  int start = 0;
  int end = ntracks;
  bool has_it = false;
  if (itr >= 0) {
    start = itr;
    end = itr + 1;
  }
  for (itr = start; itr < end; ++itr) {
    if ((fParticleV[itr] == trk) && (fEventV[itr] == evt) &&
        ((fNstepsV[itr] >= stp) && (fNstepsV[itr] < stp + nsteps))) {
      has_it = true;
#ifndef VECCORE_CUDA
      PrintTrack(itr, msg);
#else
      (void)msg;
#endif
      break;
    }
  }
  if (!has_it)
    return false;
  // Put breakpoint at line below
  return true;
}

#define IsNan(x)  ( ! ( x > 0 || x <= 0.0 ) )
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTrack_v::CheckTrack(int itr, const char *msg, double epsilon ) const
{
   // Ensure that values are 'sensible' - else print msg and track
   if( epsilon <= 0.0 || epsilon > 0.01 ) { epsilon = 1.e-6; }

   double x= fXposV[itr], y = fYposV[itr], z= fZposV[itr];
   bool badPosition = IsNan(x) || IsNan(y) || IsNan(z);
   const double maxRadius = 10000.0;   // Should be a property of the geometry
   const double maxRadXY  =  5000.0;   // Should be a property of the geometry

   // const double maxUnitDev =  1.0e-4;  // Deviation from unit of the norm of the direction
   double radiusXy2= x*x + y*y;
   double radius2=  radiusXy2 + z*z;
   badPosition = badPosition
                || (radiusXy2 > maxRadXY  * maxRadXY )
                || (radius2   > maxRadius * maxRadius ) ;

   const double maxUnitDev = epsilon;  // Use epsilon for max deviation of direction norm from 1.0

   double dx= fYdirV[itr], dy= fYdirV[itr], dz= fZdirV[itr];
   double dirNorm2 = dx*dx + dy*dy + dz*dz;
   bool   badDirection = std::fabs( dirNorm2 - 1.0 ) > maxUnitDev;
   if( badPosition || badDirection ) {
      static const char* errMsg[4]= { " All ok - No error. ",
                                      " Bad position.",                  // [1]
                                      " Bad direction.",                 // [2]
                                      " Bad direction and position. " }; // [3]
      int iM=0;
      if( badPosition ) { iM++; }
      if( badDirection  ) { iM += 2; }
      // if( badDirection ) {
      //   Printf( " Norm^2 direction= %f ,  Norm -1 = %g", dirNorm2, sqrt(dirNorm2)-1.0 );
      // }
      Printf("ERROR> Problem with track %d . Issue: %s. Info message: %s -- Mag^2(dir)= %9.6f Norm-1= %g",
             itr, errMsg[iM], msg, dirNorm2, sqrt(dirNorm2)-1.0 );
      PrintTrack(itr, msg);
   }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool GeantTrack_v::CheckDirection(int itr, double epsilon ) const
{
   if( epsilon <= 0.0 || epsilon > 0.001 ) { epsilon = 1.e-6; }

   double xdir= fXdirV[itr];
   double ydir= fYdirV[itr];
   double zdir= fZdirV[itr];
   double norm2= xdir * xdir + ydir * ydir + zdir * zdir;

   return ( std::fabs( norm2 - 1.0 ) <= epsilon );
}

} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA_ENABLED
#ifndef VECCORE_CUDA
          
//______________________________________________________________________________
bool ToDevice(vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> dest, cxx::GeantTrack_v *source, cudaStream_t stream) {
  // Since fPathV and fNextpathV are internal pointer, we need to fix them up.
  // assert(vecgeom::cuda::NavigationState::SizeOfInstance(fMaxDepth)
  //       == vecgeom::cxx::NavigationState::SizeOfInstance(fMaxDepth) );

  size_t bufferOffset = GeantTrack::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf());
  long offset = ((const char *)dest.GetPtr() + bufferOffset) - (const char *)source->Buffer();
  for (int hostIdx = 0; hostIdx < source->GetNtracks(); ++hostIdx) {
    // Technically this offset is a 'guess' and depends on the
    // host (cxx) and device (cuda) GeantTrack_v to be strictly aligned.
    if (source->fPathV[hostIdx])
      source->fPathV[hostIdx] = (VolumePath_t *)(((char *)source->fPathV[hostIdx]) + offset);
    if (source->fNextpathV[hostIdx])
      source->fNextpathV[hostIdx] = (VolumePath_t *)(((char *)source->fNextpathV[hostIdx]) + offset);
  }
  // const char* destBuf =  ((const char*)dest.GetPtr() + bufferOffset;
  // const char* sourBuf =  (const char*)source->Buffer();
  // for(int hostIdx = 0; hostIdx < source->GetNtracks(); ++hostIdx ) {
  //    fprintf(stderr,"Track[%d] : val=%p diff=%p off=%p\n", hostIdx, source->fPathV[hostIdx],
  //            ((const char*)source->fPathV[hostIdx]) - destBuf,  ((const char*)source->fPathV[hostIdx]) - offset);
  // }

  assert(((void *)source) == ((void *)(&(source->fNtracks))));

  // fprintf(stderr,"Posting the copy from host=%p to device=%p and size=%ld\n",
  //        source->Buffer(),
  //        ((char*)dest.GetPtr()) + bufferOffset,
  //        source->BufferSize());
  // fMaxtracks, fMaxDepth and fBufSize ought to be invariant.
  GEANT_CUDA_ERROR(cudaMemcpyAsync(((char*)dest.GetPtr()) + bufferOffset,
                                   source->Buffer(),
                                   source->BufferSize(),
                                   cudaMemcpyHostToDevice, stream));
  // Copy stream->fInputBasket->fNtracks, stream->fInputBasket->fNselected, stream->fInputBasket->fCompact, stream->fInputBasket->fMixed
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest,
                                   source,
                                   sizeof(int)*2+sizeof(bool)*2,
                                   cudaMemcpyHostToDevice, stream));

  return true;
}
          
//______________________________________________________________________________
void FromDeviceConversion(cxx::GeantTrack_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> source) {
  size_t bufferOffset = GeantTrack::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf());
  // Since fPathV and fNextpathV are internal pointer, we need to fix them up.
  // assert(vecgeom::cuda::NavigationState::SizeOfInstance(fMaxDepth)
  //        == vecgeom::cxx::NavigationState::SizeOfInstance(fMaxDepth) );

  long offset = ((const char *)dest->Buffer()) - (((const char *)source.GetPtr()) + bufferOffset);
  for (int hostIdx = 0; hostIdx < dest->GetNtracks(); ++hostIdx) {
    // Technically this offset is a 'guess' and depends on the
    // host (cxx) and device (cuda) GeantTrack_v to be strictly aligned.
    if (dest->fPathV[hostIdx])
      dest->fPathV[hostIdx] = (VolumePath_t *)(((char *)dest->fPathV[hostIdx]) + offset);
    if (dest->fNextpathV[hostIdx])
      dest->fNextpathV[hostIdx] = (VolumePath_t *)(((char *)dest->fNextpathV[hostIdx]) + offset);
  }
}

//______________________________________________________________________________          
bool FromDevice(cxx::GeantTrack_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> source, cudaStream_t stream) {
  size_t bufferOffset = GeantTrack::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrack_v>::SizeOf());
  // fMaxtracks, fMaxDepth and fBufSize ought to be invariant.
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest,
                                   source.GetPtr(),
                                   sizeof(int)*2+sizeof(bool)*2,
                                   cudaMemcpyDeviceToHost, stream));
  // fprintf(stderr,"Posting the copy from device=%p to host=%p and size=%lu\n",
  //        ((char*)source.GetPtr()) + bufferOffset,
  //        dest->Buffer(),
  //        dest->BufferSize());
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest->Buffer(),
                                   ((char*)source.GetPtr()) + bufferOffset,
                                   dest->BufferSize(),
                                   cudaMemcpyDeviceToHost, stream));
  return true;
}
#endif
#endif

} // Geant
