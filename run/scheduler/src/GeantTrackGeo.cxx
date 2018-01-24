#include "GeantTrackGeo.h"

#include "globals.h"
#include "Geant/Error.h"
#include "Geant/math_wrappers.h"
#include <execinfo.h>

#include "TransportManager.h"

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

#include "GeantScheduler.h"

#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "FieldLookup.h"
// #include "ConstBzFieldHelixStepper.h"
// #include "ConstVecFieldHelixStepper.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h>
#else
#include "mm_malloc.h"
#endif
#include <cassert>

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

#ifdef USE_VECGEOM_NAVIGATOR
using namespace VECGEOM_NAMESPACE;
#endif

//______________________________________________________________________________
GeantTrackGeo_v::GeantTrackGeo_v()
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0),
      fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0) {
  // Dummy ctor.
}

//______________________________________________________________________________
GeantTrackGeo_v::GeantTrackGeo_v(int size)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0),
      fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0) {
  // Constructor with maximum capacity.
  Resize(size);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTrackGeo_v *GeantTrackGeo_v::MakeInstanceAt(void *addr, unsigned int nTracks) {
  return new (addr) GeantTrackGeo_v(addr, nTracks);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTrackGeo_v::GeantTrackGeo_v(void *addr, unsigned int nTracks)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0),
      fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0) {

  // Constructor with maximum capacity.
  fBuf = ((char *)addr) + RoundUpAlign(sizeof(GeantTrackGeo_v));
  fBufSize = BufferSize(nTracks);
  memset(fBuf, 0, fBufSize);
  AssignInBuffer(fBuf, nTracks);
}

//______________________________________________________________________________
GeantTrackGeo_v::~GeantTrackGeo_v() {
  // Destructor.
  _mm_free(fBuf);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
void GeantTrackGeo_v::AssignInBuffer(char *buff, int size) {
  // Assign all internal class arrays in the supplied buffer, padded by supplied
  // size.

  const int size_doublen = size * sizeof(double);
  char *buf = buff;
  fOriginalV = (GeantTrack **)buf;
  buf += size * sizeof(GeantTrack *);
  fIdV = (size_t *)buf;
  buf += size * sizeof(size_t);
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
  fPstepV = (double *)buf;
  buf += size_doublen;
  fStepV = (double *)buf;
  buf += size_doublen;
  fSnextV = (double *)buf;
  buf += size_doublen;
  fSafetyV = (double *)buf;
  buf += size_doublen;
  fCompSafetyV = (bool *)buf;
//  buf += size_booln;
}

//______________________________________________________________________________
void GeantTrackGeo_v::CopyToBuffer(char *buff, int size) {
  // Copy existing track arrays into new buffer, padded by supplied size
  const int size_double = fNtracks * sizeof(double);
  const int size_doublen = size * sizeof(double);
  char *buf = buff;
  memcpy(buf, fOriginalV, size*sizeof(GeantTrack*));
  fOriginalV = (GeantTrack **)buf;
  buf += size*sizeof(GeantTrack*);
  memcpy(buf, fIdV, size*sizeof(size_t));
  fIdV = (size_t *)buf;
  buf += size*sizeof(size_t);
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
  memcpy(buf, fCompSafetyV, fNtracks * sizeof(bool));
  fCompSafetyV = (bool *)buf;
//  buf += size * sizeof(bool);
}

//______________________________________________________________________________
bool GeantTrackGeo_v::IsNormalized(int itr, double tolerance) const {
  // Check if track direction is normalized within tolerance
  double norm = fXdirV[itr] * fXdirV[itr] + fYdirV[itr] * fYdirV[itr] + fZdirV[itr] * fZdirV[itr];
  if (fabs(1. - norm) > tolerance)
    return false;
  return true;
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t GeantTrackGeo_v::BufferSize(size_t nTracks) {
  // return the contiguous memory size needed to hold a GeantTrackGeo's data
  size_t size = RoundUpAlign(nTracks);
  return size * sizeof(GeantTrackGeo);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t GeantTrackGeo_v::SizeOfInstance(size_t nTracks) {
  // return the contiguous memory size needed to hold a GeantTrackGeo

  return RoundUpAlign(sizeof(GeantTrackGeo_v))+BufferSize(nTracks);
}

//______________________________________________________________________________
void GeantTrackGeo_v::Resize(int newsize) {
  // Resize the container.
  int size = RoundUpAlign(newsize);
  if (size < GetNtracks()) {
    Geant::Error("Resize","%s","Cannot resize to less than current track content");
    return;
  }
  fBufSize = BufferSize(size);

  char *buf = (char *)_mm_malloc(fBufSize, GEANT_ALIGN_PADDING);
  memset(buf, 0, fBufSize);
  fMaxtracks = size;
  if (!fBuf) {
    // All arrays are contiguous in a single buffer and aligned with the
    // same padding GEANT_ALIGN_PADDING
    fBuf = buf;
    AssignInBuffer(buf, size);
  } else {
    // Resize container
    CopyToBuffer(buf, size);
    _mm_free(fBuf);
    fBuf = buf;
  }
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
int GeantTrackGeo_v::AddTracks(TrackVec_t const &array) {
  // Add all tracks from a vector into the SOA array. 
  // Returns the number of tracks after the operation.

  if ((fNtracks + (int)array.size()) >= fMaxtracks) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    Resize(Math::Max<int>(2 * fMaxtracks, fNtracks + array.size()));
#else
    printf("Error in GeantTrackGeo::AddTrack, resizing is not supported in device code\n");
#endif
  }

  for (auto track : array) {
    fOriginalV[fNtracks] = track;
    fIdV[fNtracks] = fNtracks;
    fXposV[fNtracks] = track->X();
    fYposV[fNtracks] = track->Y();
    fZposV[fNtracks] = track->Z();
    fXdirV[fNtracks] = track->Dx();
    fYdirV[fNtracks] = track->Dy();
    fZdirV[fNtracks] = track->Dz();
    fPstepV[fNtracks] = track->GetPstep();
    fStepV[fNtracks] = track->GetStep();
    fSnextV[fNtracks] = track->GetSnext();
    fSafetyV[fNtracks] = track->GetSafety();
    fCompSafetyV[fNtracks] = !track->Boundary();
    fNtracks++;
  }
  return fNtracks;  
}

} // GEANT_IMPL_NAMESPACE
} // Geant
