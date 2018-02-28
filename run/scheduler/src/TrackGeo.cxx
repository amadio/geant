#include "TrackGeo.h"

#include "Geant/Error.h"
#include "Geant/math_wrappers.h"
#include <execinfo.h>

#include "TransportManager.h"

#include "Geant/ScalarNavInterfaceVG.h"
#include "Geant/ScalarNavInterfaceVGM.h"
#include "Geant/VectorNavInterface.h"
#include "navigation/VNavigator.h"
#include "navigation/SimpleNavigator.h"
#include "navigation/ABBoxNavigator.h"
#include "volumes/PlacedVolume.h" // equivalent of TGeoNode
#include "base/Vector3D.h"
#include "base/Transformation3D.h"
#include "base/Global.h"
#include "management/GeoManager.h"
#include "base/SOA3D.h"

#include "Geant/WorkloadManager.h"

#include "Geant/TaskData.h"

#include "Geant/GUFieldPropagatorPool.h"
#include "Geant/GUFieldPropagator.h"
#include "Geant/FieldLookup.h"
// #include "Geant/ConstBzFieldHelixStepper.h"
// #include "ConstVecFieldHelixStepper.h"

#ifdef __INTEL_COMPILER
#include <immintrin.h>
#else
#include "mm_malloc.h"
#endif
#include <cassert>

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace VECGEOM_NAMESPACE;

//______________________________________________________________________________
TrackGeo_v::TrackGeo_v()
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0), fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0),
      fXdirV(0), fYdirV(0), fZdirV(0), fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0)
{
  // Dummy ctor.
}

//______________________________________________________________________________
TrackGeo_v::TrackGeo_v(int size)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0), fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0),
      fXdirV(0), fYdirV(0), fZdirV(0), fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0)
{
  // Constructor with maximum capacity.
  Resize(size);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
TrackGeo_v *TrackGeo_v::MakeInstanceAt(void *addr, unsigned int nTracks)
{
  return new (addr) TrackGeo_v(addr, nTracks);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
TrackGeo_v::TrackGeo_v(void *addr, unsigned int nTracks)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fBuf(0), fOriginalV(0), fIdV(0), fXposV(0), fYposV(0), fZposV(0),
      fXdirV(0), fYdirV(0), fZdirV(0), fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fCompSafetyV(0)
{

  // Constructor with maximum capacity.
  fBuf     = ((char *)addr) + RoundUpAlign(sizeof(TrackGeo_v));
  fBufSize = BufferSize(nTracks);
  memset(fBuf, 0, fBufSize);
  AssignInBuffer(fBuf, nTracks);
}

//______________________________________________________________________________
TrackGeo_v::~TrackGeo_v()
{
  // Destructor.
  _mm_free(fBuf);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
void TrackGeo_v::AssignInBuffer(char *buff, int size)
{
  // Assign all internal class arrays in the supplied buffer, padded by supplied
  // size.

  const int size_doublen = size * sizeof(double);
  char *buf              = buff;
  fOriginalV             = (Track **)buf;
  buf += size * sizeof(Track *);
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
void TrackGeo_v::CopyToBuffer(char *buff, int size)
{
  // Copy existing track arrays into new buffer, padded by supplied size
  const int size_double  = fNtracks * sizeof(double);
  const int size_doublen = size * sizeof(double);
  char *buf              = buff;
  memcpy(buf, fOriginalV, size * sizeof(Track *));
  fOriginalV = (Track **)buf;
  buf += size * sizeof(Track *);
  memcpy(buf, fIdV, size * sizeof(size_t));
  fIdV = (size_t *)buf;
  buf += size * sizeof(size_t);
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
bool TrackGeo_v::IsNormalized(int itr, double tolerance) const
{
  // Check if track direction is normalized within tolerance
  double norm = fXdirV[itr] * fXdirV[itr] + fYdirV[itr] * fYdirV[itr] + fZdirV[itr] * fZdirV[itr];
  if (fabs(1. - norm) > tolerance) return false;
  return true;
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t TrackGeo_v::BufferSize(size_t nTracks)
{
  // return the contiguous memory size needed to hold a TrackGeo's data
  size_t size = RoundUpAlign(nTracks);
  return size * sizeof(TrackGeo);
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t TrackGeo_v::SizeOfInstance(size_t nTracks)
{
  // return the contiguous memory size needed to hold a TrackGeo

  return RoundUpAlign(sizeof(TrackGeo_v)) + BufferSize(nTracks);
}

//______________________________________________________________________________
void TrackGeo_v::Resize(int newsize)
{
  // Resize the container.
  int size = RoundUpAlign(newsize);
  if (size < GetNtracks()) {
    geant::Error("Resize", "%s", "Cannot resize to less than current track content");
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
int TrackGeo_v::AddTracks(TrackVec_t const &array)
{
  // Add all tracks from a vector into the SOA array.
  // Returns the number of tracks after the operation.

  if ((fNtracks + (int)array.size()) >= fMaxtracks) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    Resize(Math::Max<int>(2 * fMaxtracks, fNtracks + array.size()));
#else
    printf("Error in TrackGeo::AddTrack, resizing is not supported in device code\n");
#endif
  }

  for (auto track : array) {
    fOriginalV[fNtracks]   = track;
    fIdV[fNtracks]         = fNtracks;
    fXposV[fNtracks]       = track->X();
    fYposV[fNtracks]       = track->Y();
    fZposV[fNtracks]       = track->Z();
    fXdirV[fNtracks]       = track->Dx();
    fYdirV[fNtracks]       = track->Dy();
    fZdirV[fNtracks]       = track->Dz();
    fPstepV[fNtracks]      = track->GetPstep();
    fStepV[fNtracks]       = track->GetStep();
    fSnextV[fNtracks]      = track->GetSnext();
    fSafetyV[fNtracks]     = track->GetSafety();
    fCompSafetyV[fNtracks] = !track->Boundary();
    fNtracks++;
  }
  return fNtracks;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
