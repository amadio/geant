#include "GeantTrackGeo.h"

#include "globals.h"
#include "Geant/Error.h"
#include <execinfo.h>

#ifdef USE_VECGEOM_NAVIGATOR
#pragma message("Compiling against VecGeom")
#include "ScalarNavInterfaceVG.h"
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
#include "backend/Backend.h"
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
#pragma message("Compiling against TGeo")
#include "ScalarNavInterfaceTGeo.h"
#include <iostream>
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#endif

#include "WorkloadManager.h"

#include "GeantTaskData.h"
#include "ConstFieldHelixStepper.h"
#include "GeantScheduler.h"

// #ifdef  RUNGE_KUTTA
#pragma message("Compiling using Runge-Kutta for integration")
#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
// #endif

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
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fVPstart(0), fBuf(0),
      fOriginalV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fBoundaryV(0), fPathV(0), fNextpathV(0) {
  // Dummy ctor.
}

//______________________________________________________________________________
GeantTrackGeo_v::GeantTrackGeo_v(int size)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fVPstart(0), fBuf(0),
      fOriginalV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fBoundaryV(0), fPathV(0), fNextpathV(0) {
  // Constructor with maximum capacity.
  Resize(size);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrackGeo_v *GeantTrackGeo_v::MakeInstanceAt(void *addr, unsigned int nTracks) {
  return new (addr) GeantTrackGeo(addr, nTracks);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
GeantTrackGeo_v::GeantTrackGeo_v(void *addr, unsigned int nTracks)
    : fNtracks(0), fMaxtracks(0), fBufSize(0), fVPstart(0), fBuf(0),
      fOriginalV(0), fXposV(0), fYposV(0), fZposV(0), fXdirV(0), fYdirV(0), fZdirV(0),
      fPstepV(0), fStepV(0), fSnextV(0), fSafetyV(0), fBoundaryV(0), fPathV(0), fNextpathV(0) {

  // Constructor with maximum capacity.
  fBuf = ((char *)addr) + round_up_align(sizeof(GeantTrackGeo_v));
  fBufSize = BufferSize(nTracks);
  memset(fBuf, 0, fBufSize);
  AssignInBuffer(fBuf, nTracks);
  memset(fPathV, 0, nTracks * sizeof(VolumePath_t *));
  memset(fNextpathV, 0, nTracks * sizeof(VolumePath_t *));
}

//______________________________________________________________________________
GeantTrackGeo_v::~GeantTrackGeo_v() {
  // Destructor.
  for (auto i = 0; i < fNtracks; ++i) {
    VolumePath_t::ReleaseInstance(fPathV[i]);
    VolumePath_t::ReleaseInstance(fNextpathV[i]);
  }
  _mm_free(fBuf);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo_v::AssignInBuffer(char *buff, int size) {
  // Assign all internal class arrays in the supplied buffer, padded by supplied
  // size.

  const int size_intn = size * sizeof(int);
  const int size_doublen = size * sizeof(double);
  const int size_booln = size * sizeof(bool);
  char *buf = buff;
  fOriginalV = (GeantTrack *)buf;
  buf += size * sizeof(GeantTrack *);
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
  fBoundaryV = (bool *)buf;
  buf += size_booln;
  fPathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  fNextpathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
}

//______________________________________________________________________________
void GeantTrackGeo_v::CopyToBuffer(char *buff, int size) {
  // Copy existing track arrays into new buffer, padded by supplied size
  const int size_int = fNtracks * sizeof(int);
  const int size_double = fNtracks * sizeof(double);
  const int size_intn = size * sizeof(int);
  const int size_doublen = size * sizeof(double);
  char *buf = buff;
  memcpy(buf, fOriginalV, size*sizeof(GeantTrack*));
  fOriginalV = (GeantTrack **)buf;
  buf += size*sizeof(GeantTrack*);
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
  memcpy(buf, fBoundaryV, ntracks * sizeof(bool));
  fBoundaryV = (bool *)buf;
  buf += size * sizeof(bool);
  memcpy(buf, fPathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **pathV = (VolumePath_t **)buf;
  buf += size * sizeof(VolumePath_t *);
  memcpy(buf, fNextpathV, ntracks*sizeof(VolumePath_t*));
  VolumePath_t **nextpathV = (VolumePath_t **)buf;
  //buf += size * sizeof(VolumePath_t *);
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
GEANT_CUDA_BOTH_CODE
size_t GeantTrackGeo_v::BufferSize(size_t nTracks) {
  // return the contiguous memory size needed to hold a GeantTrackGeo's data
  size_t size = round_up_align(nTracks);
  return size * sizeof(GeantTrackGeo);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
size_t GeantTrackGeo_v::SizeOfInstance(size_t nTracks) {
  // return the contiguous memory size needed to hold a GeantTrackGeo

  return round_up_align(sizeof(GeantTrackGeo_v))+BufferSize(nTracks);
}

//______________________________________________________________________________
void GeantTrackGeo_v::Resize(int newsize) {
  // Resize the container.
  int size = round_up_align(newsize);
  if (size < GetNtracks()) {
    Geant::Error("Resize","Cannot resize to less than current track content");
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
    memset(fPathV, 0, size * sizeof(VolumePath_t *));
    memset(fNextpathV, 0, size * sizeof(VolumePath_t *));
  } else {
    // Resize container
    CopyToBuffer(buf, size);
    _mm_free(fBuf);
    fBuf = buf;
  }
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo_v::AddTracks(TrackArray_t const &array) {
  // Add all tracks from a vector into the SOA array. 
  // Returns the number of tracks after the operation.

  if (fNtracks + array.size() >= fMaxtracks) {
#ifndef GEANT_CUDA_DEVICE_BUILD
    Resize(Max(2 * fMaxtracks, fNtracks + array.size()));
#else
    printf("Error in GeantTrackGeo::AddTrack, resizing is not supported in device code\n");
#endif
  }

  for (auto track : array) {
    fOriginalV[fNtracks] = track;
    fXposV[fNtracks] = track->fXpos;
    fYposV[fNtracks] = track->fYpos;
    fZposV[fNtracks] = track->fZpos;
    fXdirV[fNtracks] = track->fXdir;
    fYdirV[fNtracks] = track->fYdir;
    fZdirV[fNtracks] = track->fZdir;
    fPstepV[fNtracks] = track->fPstep;
    fStepV[fNtracks] = track->fStep;
    fSnextV[fNtracks] = track->fSnext;
    fSafetyV[fNtracks] = track->fSafety;
    fBoundaryV[fNtracks] = track->fBoundary;
    fPathV[fNtracks] = track->fPath;
    fNextpathV[fNtracks] = track->fNextpath;
    fNtracks++;
  }
  return fNtracks;  
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo_v::UpdateOriginalTrack(int itr) const {
  // Update the original track itr.
  GeantTrackGeo &track = *fOriginals[i];
  track.fXpos = fXposV[i];
  track.fYpos = fYposV[i];
  track.fZpos = fZposV[i];
  track.fXdir = fXdirV[i];
  track.fYdir = fYdirV[i];
  track.fZdir = fZdirV[i];
  track.fPstep = fPstepV[i];
  track.fStep = fStepV[i];
  track.fSnext = fSnextV[i];
  track.fSafety = fSafetyV[i];
  track.fBoundary = fBoundaryV[i];
  // The path and nextpath members are already pointing to the storage used by the original track
}

//______________________________________________________________________________
int GeantTrackGeo::PropagateStraight(int ntracks, double *crtstep) {
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
void GeantTrackGeo::PropagateInVolume(int ntracks, const double *crtstep, GeantTaskData *td) {
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
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::PropagateInVolumeSingle(int i, double crtstep, GeantTaskData * td) {
  // Propagate the selected track with crtstep value. The method is to be called
  // only with  charged tracks in magnetic field.The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)

   // Double_t c = 0.;
   // const Double_t *point = 0;
   // const Double_t *newdir = 0;

   bool useRungeKutta;
#ifdef GEANT_CUDA_DEVICE_BUILD
   const double bmag = gPropagator_fBmag;
   constexpr auto gPropagator_fUseRK = false; // Temporary work-around until actual implementation ..
   useRungeKutta= gPropagator_fUseRK;   //  Something like this is needed - TBD
#else
   const double bmag = gPropagator->fBmag;
   useRungeKutta= gPropagator->fUseRungeKutta;
#endif

   // static unsigned long icount= 0;
   // if( icount++ < 2 )  std::cout << " PropagateInVolumeSingle: useRungeKutta= " << useRungeKutta << std::endl;

// #ifdef RUNGE_KUTTA
#ifndef GEANT_CUDA_DEVICE_BUILD
   GUFieldPropagator *fieldPropagator = nullptr;
   if( useRungeKutta ){
      // Initialize for the current thread -- move to GeantPropagator::Initialize()
      static GUFieldPropagatorPool* fieldPropPool= GUFieldPropagatorPool::Instance();
      assert( fieldPropPool );

      fieldPropagator = fieldPropPool->GetPropagator(td->fTid);
      assert( fieldPropagator );
   }
#endif

  // Reset relevant variables
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
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
// alternative code with lean stepper would be:
// ( stepper header has to be included )

  using ThreeVector = vecgeom::Vector3D<double>;
  // typedef vecgeom::Vector3D<double>  ThreeVector;   
  ThreeVector Position(fXposV[i], fYposV[i], fZposV[i]);
  ThreeVector Direction(fXdirV[i], fYdirV[i], fZdirV[i]);
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

  if( useRungeKutta ) {
#ifndef GEANT_NVCC
     fieldPropagator->DoStep(Position,    Direction,    fChargeV[i], fPV[i], crtstep,
                             PositionNew, DirectionNew);
#endif
  } else {
     // Old - constant field
     Geant::ConstBzFieldHelixStepper stepper(bmag);
     stepper.DoStep<ThreeVector,double,int>(Position,    Direction,    fChargeV[i], fPV[i], crtstep,
                                         PositionNew, DirectionNew);
  }

  fXposV[i] = PositionNew.x();
  fYposV[i] = PositionNew.y();
  fZposV[i] = PositionNew.z();

  //  maybe normalize direction here  // Math::Normalize(dirnew);
  DirectionNew = DirectionNew.Unit();   
  fXdirV[i] = DirectionNew.x();
  fYdirV[i] = DirectionNew.y();
  fZdirV[i] = DirectionNew.z();

#if 0
  ThreeVector SimplePosition = Position + crtstep * Direction;
  // double diffpos2 = (PositionNew - Position).Mag2();
  double diffpos2 = (PositionNew - SimplePosition).Mag2();
  //   -- if (Math::Sqrt(diffpos)>0.01*crtstep) {     
  const double drift= 0.01*crtstep;
  if ( diffpos2>drift*drift ){
      double diffpos= Math::Sqrt(diffpos2);
      // Geant::Print("PropagateInVolumeSingle","relative difference in pos = %g", diffpos/crtstep);
      Geant::Print("PropagateInVolumeSingle","difference in pos = %g (abs) %g (relative) , step= %g",
                   diffpos, diffpos/crtstep, crtstep);
  }
#endif
}

#ifdef USE_VECGEOM_NAVIGATOR
void GeantTrackGeo::CheckLocationPathConsistency(int itr) const {
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
int GeantTrackGeo::SortByStatus(TrackStatus_t status) {
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
int GeantTrackGeo::SortByLimitingDiscreteProcess() {
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
int GeantTrackGeo::RemoveByStatus(TrackStatus_t status, GeantTrackGeo &output) {
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
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::PrintTrack(int itr, const char *msg) const {
  // Print info for a given track
  const char *status[8] = {"alive", "killed", "inflight", "boundary", "exitSetup", "physics", "postponed", "new"};
#ifdef USE_VECGEOM_NAVIGATOR
  Geant::Print(msg,
      "== Track %d: evt=%d slt=%d part=%d pdg=%d gVc=%d chg=%d proc=%d nstp=%d spc=%d status=%s mass=%g "
      "xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g time=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n",
      itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fPDGV[itr], fGVcodeV[itr],
      fChargeV[itr], fProcessV[itr], fNstepsV[itr], (int)fSpeciesV[itr], status[int(fStatusV[itr])],
      fMassV[itr], fXposV[itr], fYposV[itr], fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr],
      fTimeV[itr], fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fNintLenV[itr], fIntLenV[itr], fBoundaryV[itr]);
  
#ifndef GEANT_NVCC
  fPathV[itr]->Print();
  fNextpathV[itr]->Print();
#endif
#else
  TString path;
  fPathV[itr]->GetPath(path);
  TString nextpath;
  fNextpathV[itr]->GetPath(nextpath);

  Geant::Print(msg, "== Track %d: evt=%d slt=%d part=%d pdg=%d gVc=%d eind=%d chg=%d proc=%d nstp=%d "
         "spc=%d status=%s mass=%g xpos=%g ypos=%g zpos=%g xdir=%g ydir=%g zdir=%g mom=%g ene=%g "
         "time=%g edep=%g pstp=%g stp=%g snxt=%g saf=%g nil=%g ile=%g bdr=%d\n pth=%s npth=%s\n",
         itr, fEventV[itr], fEvslotV[itr], fParticleV[itr], fPDGV[itr], fEindexV[itr], fGVcodeV[itr],
         fChargeV[itr], fProcessV[itr], fNstepsV[itr], (int)fSpeciesV[itr], status[int(fStatusV[itr])],
         fMassV[itr], fXposV[itr], fYposV[itr], fZposV[itr], fXdirV[itr], fYdirV[itr], fZdirV[itr], fPV[itr], fEV[itr],
         fTimeV[itr], fEdepV[itr], fPstepV[itr], fStepV[itr], fSnextV[itr], fSafetyV[itr], fNintLenV[itr], fIntLenV[itr], fBoundaryV[itr], path.Data(),
         nextpath.Data());
#endif
}

//______________________________________________________________________________
void GeantTrackGeo::PrintTracks(const char *msg) const {
  // Print all tracks
  int ntracks = GetNtracks();
  Geant::Print(msg,"");
  for (int i = 0; i < ntracks; i++)
    PrintTrack(i);
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::ComputeTransportLength(int ntracks, GeantTaskData *td) {
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
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::ComputeTransportLengthSingle(int itr, GeantTaskData *td) {
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
TransportAction_t GeantTrackGeo::PostponedAction(int ntracks) const {
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
int GeantTrackGeo::PropagateTracks(GeantTaskData *td) {
  // Propagate the ntracks in the current volume with their physics steps (already
  // computed)
  // Vectors are pushed downstream when efficient.
  GeantTrackGeo &output = *td->fTransported;
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
  GeantPropagator *prop = GeantPropagator::Instance();
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "PropagateTracks");
#endif
  ComputeTransportLength(ntracks, td);
//         Printf("====== After ComputeTransportLength:");
//         PrintTracks();
#ifdef BUG_HUNT
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterCompTransLen");
#endif

  int itr = 0;
  int icrossed = 0;
  int nsel = 0;
  double lmax;
  const double eps = 1.E-2; // 100 micron
  const double bmag = gPropagator->fBmag;

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
    if (fChargeV[itr] == 0 || bmag < 1.E-10) {
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
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateTracks");
#endif
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo::PropagateSingleTrack(int itr, GeantTaskData *td, int stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  int icrossed = 0;
  double step, lmax;
  const double eps = 1.E-2; // 1 micron
#ifdef GEANT_CUDA_DEVICE_BUILD
  const double bmag = gPropagator_fBmag;
#else
  const double bmag = gPropagator->fBmag;
#endif
// Compute transport length in geometry, limited by the physics step
#ifdef BUG_HUNT
  GeantPropagator *prop = GeantPropagator::Instance();
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "PropagateSingle", itr);
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
    if (fChargeV[itr] == 0 || bmag < 1.E-10) {
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
      BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateSingleNeutral",
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
    step = (fBoundaryV[itr]) ? Math::Min<double>(lmax, Math::Max<double>(fSnextV[itr], 1.E-4)) : Math::Min<double>(lmax, fPstepV[itr]);
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
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateSingle", itr);
#endif
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo::PropagateTracksScalar(GeantTaskData *td, int stage) {
  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

#ifndef GEANT_CUDA_DEVICE_BUILD
  GeantTrackGeo &output = *td->fTransported;
#endif
  int icrossed = 0;
  int ntracks = GetNtracks();
  for (int itr = 0; itr < ntracks; itr++) {
    icrossed += PropagateSingleTrack(itr, td, stage);
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
double GeantTrackGeo::Curvature(int i) const {
  // Curvature assuming constant field is along Z
  constexpr double kB2C = -0.299792458e-3;
  constexpr double tiny = 1.E-30;
#ifdef GEANT_CUDA_DEVICE_BUILD
  const double bmag = gPropagator_fBmag;
#else
  const double bmag = gPropagator->fBmag;
#endif
  return fabs(kB2C * fChargeV[i] * bmag / (Pt(i) + tiny));
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
double GeantTrackGeo::SafeLength(int i, double eps) {
  // Returns the propagation length in field such that the propagated point is
  // shifted less than eps with respect to the linear propagation.
  double c = Curvature(i);
  if (c < 1.E-10)
    return 1.E20;
  return 2. * sqrt(eps / c);
}

//______________________________________________________________________________
int GeantTrackGeo::PostponeTracks(GeantTrackGeo &output) {
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
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo::PostponeTrack(int itr, GeantTrackGeo &output) {
  // Postpone transport of a track and copy it to the output.
  // Returns where in the output the track was added.

  fStatusV[itr] = kPostponed;
  // Move these tracks to the output container
  int new_itr = output.AddTrack(*this, itr, true);
  MarkRemoved(itr);
  return new_itr;
}


//______________________________________________________________________________
Volume_t const*GeantTrackGeo::GetNextVolume(int i) const {
  // Next volume the track is getting into
#ifdef USE_VECGEOM_NAVIGATOR
  return fNextpathV[i]->Top()->GetLogicalVolume();
#else
  return fNextpathV[i]->GetCurrentNode()->GetVolume();
#endif
}

//______________________________________________________________________________
Volume_t const*GeantTrackGeo::GetVolume(int i) const {
  // Current volume the track is into
#ifdef USE_VECGEOM_NAVIGATOR
  return fPathV[i]->Top()->GetLogicalVolume();
#else
  return (fPathV[i]->GetCurrentNode()->GetVolume());
#endif
}

//______________________________________________________________________________
Material_t *GeantTrackGeo::GetMaterial(int i) const {
  // Current material the track is into
#ifdef USE_VECGEOM_NAVIGATOR
  Medium_t *med = (Medium_t *)GetVolume(i)->GetTrackingMediumPtr();
#else
  Medium_t *med = (Medium_t *)GetVolume(i)->GetMedium();
#endif
  // TODO: better to use assert
  if (!med)
    return nullptr;
  return med->GetMaterial();
}

//______________________________________________________________________________
#ifdef USE_VECGEOM_NAVIGATOR
bool GeantTrackGeo::CheckNavConsistency(int /*itr*/) {
  // TO IMPLEMENT WIRH VECGEOM
#else
bool GeantTrackGeo::CheckNavConsistency(int itr) {
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
GEANT_CUDA_BOTH_CODE
bool GeantTrackGeo::BreakOnStep(int evt, int trk, int stp, int nsteps, const char *msg, int itr) {
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
#ifndef GEANT_NVCC
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

} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA
#ifndef GEANT_NVCC

bool ToDevice(vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo> dest, cxx::GeantTrackGeo *source, cudaStream_t stream) {
  // Since fPathV and fNextpathV are internal pointer, we need to fix them up.
  // assert(vecgeom::cuda::NavigationState::SizeOfInstance(fMaxDepth)
  //       == vecgeom::cxx::NavigationState::SizeOfInstance(fMaxDepth) );

  size_t bufferOffset = GeantTrackGeo::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrackGeo>::SizeOf());
  long offset = ((const char *)dest.GetPtr() + bufferOffset) - (const char *)source->Buffer();
  for (int hostIdx = 0; hostIdx < source->GetNtracks(); ++hostIdx) {
    // Technically this offset is a 'guess' and depends on the
    // host (cxx) and device (cuda) GeantTrackGeo to be strictly aligned.
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

  fprintf(stderr,"Posting the copy from host=%p to device=%p and size=%ld\n",
          source->Buffer(),
          ((char*)dest.GetPtr()) + bufferOffset,
          source->BufferSize());
  // fMaxtracks, fMaxDepth and fBufSize ought to be invariant.
  GEANT_CUDA_ERROR(cudaMemcpyAsync(((char*)dest.GetPtr()) + bufferOffset,
                                   source->Buffer(),
                                   source->BufferSize(),
                                   cudaMemcpyHostToDevice, stream));
  // Copy stream->fInputBasket->fNtracks, stream->fInputBasket->fNselected, stream->fInputBasket->fCompact, stream->fInputBasket->fMixed
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest,
                                   source,
                                   sizeof(int)*2+sizeof(Bool_t)*2,
                                   cudaMemcpyHostToDevice, stream));

  return true;
}

void FromDeviceConversion(cxx::GeantTrackGeo *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo> source) {
  size_t bufferOffset = GeantTrackGeo::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrackGeo>::SizeOf());
  // Since fPathV and fNextpathV are internal pointer, we need to fix them up.
  // assert(vecgeom::cuda::NavigationState::SizeOfInstance(fMaxDepth)
  //        == vecgeom::cxx::NavigationState::SizeOfInstance(fMaxDepth) );

  long offset = ((const char *)dest->Buffer()) - (((const char *)source.GetPtr()) + bufferOffset);
  for (int hostIdx = 0; hostIdx < dest->GetNtracks(); ++hostIdx) {
    // Technically this offset is a 'guess' and depends on the
    // host (cxx) and device (cuda) GeantTrackGeo to be strictly aligned.
    if (dest->fPathV[hostIdx])
      dest->fPathV[hostIdx] = (VolumePath_t *)(((char *)dest->fPathV[hostIdx]) + offset);
    if (dest->fNextpathV[hostIdx])
      dest->fNextpathV[hostIdx] = (VolumePath_t *)(((char *)dest->fNextpathV[hostIdx]) + offset);
  }
}

bool FromDevice(cxx::GeantTrackGeo *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo> source, cudaStream_t stream) {
  size_t bufferOffset = GeantTrackGeo::round_up_align(vecgeom::cxx::DevicePtr<Geant::cuda::GeantTrackGeo>::SizeOf());
  // fMaxtracks, fMaxDepth and fBufSize ought to be invariant.
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest,
                                   source.GetPtr(),
                                   sizeof(int)*2+sizeof(Bool_t)*2,
                                   cudaMemcpyDeviceToHost, stream));
  fprintf(stderr,"Posting the copy from device=%p to host=%p and size=%lu\n",
          ((char*)source.GetPtr()) + bufferOffset,
          dest->Buffer(),
          dest->BufferSize());
  GEANT_CUDA_ERROR(cudaMemcpyAsync(dest->Buffer(),
                                   ((char*)source.GetPtr()) + bufferOffset,
                                   dest->BufferSize(),
                                   cudaMemcpyDeviceToHost, stream));
  return true;
}
#endif
#endif

} // Geant
