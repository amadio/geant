#include "TransportManager.h"

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
GEANT_CUDA_BOTH_CODE
int TransportManager::CheckSameLocation(TrackVec_t &tracks,
                                         int ntracks,
                                         GeantTaskData *td) {
// Query geometry if the location has changed for a vector of particles
// Returns number of tracks crossing the boundary.

  const int kMinVecSize = 4; // this should be retrieved from elsewhere
  TrackVec_t &output = *td->fTransported;
  int icrossed = 0;
  int nsel = ntracks;
  for (int itr = 0; itr < ntracks; itr++) {
    GeantTrack &track = *tracks[itr];
    if (track.fSafety > 1.E-10 && track.fSnext > 1.E-10) {
      // Track stays in the same volume
      track.fBoundary = false;
      nsel--;
    }
  }
  if (!nsel) return 0;
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  if (nsel < kMinVecSize) {
    // No efficient vectorization possible, do scalar treatment
#endif
    // The scalar treatment is always done in case of non-vector compilation
    for (int itr = 0; itr < tracks.size(); itr++) {
      GeantTrack &track = *tracks[itr];
      int crossed = CheckSameLocationSingle(track, td);
      if (crossed) MoveTrack(itr--, tracks, output);
      icrossed += crossed;
    }
    return icrossed;
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  }
#endif

#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  // Vector treatment
  // Copy data to SOA and dispatch for vector mode
  GeantTrackGeo_v &track_geo = td.GeoTrack();
  track_geo.Clear();
  for (int itr = 0; itr < tracks.size(); itr++) {
    GeantTrack &track = *tracks[itr];
    if (track.fSafety < 1.E-10 || track.fSnext < 1.E-10)
      track_geo.AddTrack(track);
  }
  bool *same = td->GetBoolArray(nsel);
  NavigationState *tmpstate = td->GetPath();
  VectorNavInterface::NavIsSameLocation(nsel,
                   track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
                   track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
                   (const VolumePath_t**)fPathV, fNextpathV, same, tmpstate);
  for (itr = 0; itr < nsel; itr++) {
    if (same[itr]) {
      track_geo.fBoundaryV[itr] = false;
      continue;
    }
    // Boundary crossed
    track_geo.fOriginalV[itr]->fStatus = kBoundary;
    if (track_geo.fNextpath[itr]->IsOutside())
      track_geo.fOriginalV[itr]->fStatus = kExitingSetup;
    track_geo.fBoundaryV[itr] = true;
    icrossed++;
    if (track_geo.fStepV[itr] < 1.E-8) td->fNsmall++;
  }
  track_geo.UpdateOriginalTracks();

  // Move crossing tracks to the output container
  for (int itr = 0; itr < ntracks; itr++) {
    GeantTrack &track = *tracks[itr];
    if (track.fBoundary) MoveTrack(itr--, tracks, output);
  }
  
  return icrossed;
#endif   // vector mode
}  

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int TransportManager::CheckSameLocationSingle(GeantTrack &track,
                                         GeantTaskData *td) {
// Query geometry if the location has changed for a track
// Returns number of tracks crossing the boundary (0 or 1)

  if (track.fSafety > 1.E-10 && track.fSnext > 1.E-10) {
    // Track stays in the same volume
    track.fBoundary = false;
    return 0;
  }
  
  // Track may have crossed, check it
  bool same;
#ifdef USE_VECGEOM_NAVIGATOR
  NavigationState *tmpstate = td->GetPath();
#ifdef NEW_NAVIGATION
  ScalarNavInterfaceVGM
#else
// Old navigation system
  ScalarNavInterfaceVG
#endif // NEW_NAVIGATION
    ::NavIsSameLocation(track, same, tmpstate);
#else
// ROOT navigation
  ScalarNavInterfaceTGeo
    ::NavIsSameLocation(track, same);
#endif // USE_VECGEOM_NAVIGATOR
  if (same) {
    track.fBoundary = false;
    return 0;
  }
  
  track.fBoundary = true;
  track.fStatus = kBoundary;
  if (track.fNextpath->IsOutside())
    track.fStatus = kExitingSetup;
  if (track.fStep < 1.E-8) td->fNsmall++;
  return 1;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void TransportManager::ComputeTransportLength(TrackVec_t &tracks,
                                              int ntracks,
                                              GeantTaskData *td) {
// Vector version for proposing the geometry step. All tracks have to be
// in the same volume. This may still fall back on the scalar implementation
// in case the vector size is too small.
  const int kMinVecSize = 4; // this should be retrieved from elsewhere
#ifdef USE_VECGEOM_NAVIGATOR
//#define VECTORIZED_GEOMETRY
#ifdef VECTORIZED_GEOMETRY
  // We ignore tracks for which the current safety allows for the proposed step
  // We need to count if the remaining tracks may form a vector
  int nsel = 0;
  for (auto track : tracks) 
    if (track->fSafety < track->fPstep) nsel++;
  if (nsel < kMinVecSize) {
    // Just process all the tracks in scalar mode
    for (auto track : tracks)
      ComputeTransportLengthSingle(*track, td);
    }
  } else {
    // Copy interesting tracks to SOA and process vectorized.
    // This is overhead but cannot be avoided.
    GeantTrackGeo_v &track_geo = td.GeoTrack();
    track_geo.Clear();
    for (auto track : tracks) {
      if (track->fSafety < track->fPstep)
        track_geo.AddTrack(*track);
    }
    // The vectorized SOA call
    VectorNavInterface::NavFindNextBoundaryAndStep(nsel, track_geo.fPstepV,
                 track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
                 track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
                 (const VolumePath_t **)track_geo.fPathV, track_geo.fNextpathV,
                 track_geo.fSnextV, track_geo.fSafetyV, track_geo.fBoundaryV);
  
    // Update original tracks
    track_geo->UpdateOriginalTracks();
    // Update number of calls to geometry (1 vector call + ntail scalar calls)
    td->fNsnext += 1 + ntracks%kMinVecSize;
  }
#else  // !VECTORIZED_GEOMETRY
  // Non-vectorized looped implementation of VecGeom navigation
  for (auto track : tracks) {
#ifdef NEW_NAVIGATION
    ScalarNavInterfaceVGM
#else
    ScalarNavInterfaceVG
#endif // NEW_NAVIGATION
      ::NavFindNextBoundaryAndStep(*track);
  }
  // Update number of calls to geometry (consider N scalar calls)
  td->fNsnext += ntracks;
#endif // VECTORIZED_GEOMETRY
  // perform a couple of additional checks/ set status flags and so on
  for (auto track : tracks) {
    if ((track->fNextpath->IsOutside() && track->fSnext < 1.E-6) || track->fSnext > 1.E19)
      track->fStatus = kExitingSetup;
  }
#else
  // TGeo implementation fall back on looped version
  for (auto track : tracks) {
    ComputeTransportLengthSingle(*track, td);
  }
#endif // USE_VECGEOM_NAVIGATOR 
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::ComputeTransportLengthSingle(GeantTrack &track, GeantTaskData *td) {
// Computes snext and safety for a single track. For charged tracks these are the only
// computed values, while for neutral ones the next node is checked and the boundary flag is set if
// closer than the proposed physics step.

#ifdef USE_VECGEOM_NAVIGATOR
//#define NEW_NAVIGATION
#ifdef NEW_NAVIGATION
  ScalarNavInterfaceVGM
   ::NavFindNextBoundaryAndStep(track);
#else
  ScalarNavInterfaceVG
   ::NavFindNextBoundaryAndStep(track);
#endif // NEW_NAVIGATION
#else
// ROOT geometry
  ScalarNavInterfaceTGeo
   ::NavFindNextBoundaryAndStep(track);
#endif // USE_VECGEOM_NAVIGATOR
  // Update number of calls to geometry
  td->fNsnext++;
  // if outside detector or enormous step mark particle as exiting the detector
  if (track.fNextpath->IsOutside() || track.fSnext > 1.E19)
    track.fStatus = kExitingSetup;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::PropagateInVolume(TrackVec_t &tracks,
                                      int ntracks,
                                      const double *crtstep,
                                      GeantTaskData *td) {
  // Propagate the selected tracks with crtstep values. The method is to be called
  // only with  charged tracks in magnetic field. The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)
  
  // For the moment fall back on scalar propagation
  int itr = 0;
  for (auto track : tracks) {
    PropagateInVolumeSingle(*track, crtstep[itr++], td);
  }
}


//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
void GeantTrackGeo::PropagateInVolumeSingle(GeantTrack &track, double crtstep, GeantTaskData * td) {
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
  track.fStatus = kInFlight;
  track.fPstep -= crtstep;
  if (track.fPstepV < 1.E-10) {
    track.fPstep = 0;
    track.fStatus = kPhysics;
  }
  track.fSafety -= crtstep;
  if (track.fSafety < 1.E-10)
    track.fSafety = 0;
  track.fSnext -= crtstep;
  if (track.fSnext < 1.E-10) {
    track.fSnext = 0;
    if (track.fBoundary) {
      track.fStatus = kBoundary;
    }
  }
  fStep += crtstep;
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
// alternative code with lean stepper would be:
// ( stepper header has to be included )

  using ThreeVector = vecgeom::Vector3D<double>;
  // typedef vecgeom::Vector3D<double>  ThreeVector;   
  ThreeVector Position(track.fXpos, track.fYpos, track.fZpos);
  ThreeVector Direction(track.fXdir, track.fYdir, track.fZdir);
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

  if( useRungeKutta ) {
#ifndef GEANT_NVCC
     fieldPropagator->DoStep(Position,    Direction,    track.fCharge, track.fP, crtstep,
                             PositionNew, DirectionNew);
#endif
  } else {
     // Old - constant field
     Geant::ConstBzFieldHelixStepper stepper(bmag);
     stepper.DoStep<ThreeVector,double,int>(Position,    Direction,    track.fCharge, track.fP, crtstep,
                                         PositionNew, DirectionNew);
  }

  track.fXpos = PositionNew.x();
  track.fYpos = PositionNew.y();
  track.fZpos = PositionNew.z();

  //  maybe normalize direction here  // Math::Normalize(dirnew);
  DirectionNew = DirectionNew.Unit();   
  track.fXdir = DirectionNew.x();
  track.fYdir = DirectionNew.y();
  track.fZdir = DirectionNew.z();

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

//______________________________________________________________________________
int GeantTrackGeo::PropagateTracks(TrackVec_t &tracks, GeantTaskData *td) {
  // Propagate the ntracks in the current volume with their physics steps (already
  // computed)
  // Vectors are pushed downstream when efficient.
  TrackVec_t &output = *td->fTransported;
  int ntracks = tracks.size();
  // Check if tracking the remaining tracks can be postponed
  TransportAction_t action = PostponedAction(ntracks);
  if (action == kPostpone) {
    // This needs a review. Currently kPostpone is not fired.
    PostponeTracks(tracks, output);
    // We should notify somehow the workload manager that there are postponed tracks
    return 0;
  }
  if (action != kVector)
    return PropagateTracksScalar(tracks, td);
// Compute transport length in geometry, limited by the physics step
#ifdef BUG_HUNT
  GeantPropagator *prop = GeantPropagator::Instance();
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "PropagateTracks");
#endif
  ComputeTransportLength(tracks, ntracks, td);
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
  int itr = 0;
  for (int itr=0; itr<tracks.size(); ++itr) {
    GeantTrack &track = *tracks[itr];
    // Move dead to output
    if (track.fSnext < 0) {
      Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", track.fParticle);
      track.Print();
      track.fStatus = kKilled;
    }
    if (track.fStatus == kKilled) {
      MoveTrack(itr--, tracks, output); // next track at same index after move
      continue;
    }
    // Propagate straight tracks to the precomputed location and update state,
    // ten move them to output
    if (track.fCharge == 0 || bmag < 1.E-10) {
      // Do straight propagation to physics process or boundary
      if (track.fBoundary) {
        if (track.fNextpath->IsOutside())
          track.fStatus = kExitingSetup;
        else
          track.fStatus = kBoundary;
        icrossed++;
      } else {
        track.fStatus = kPhysics;
        // Update number of steps to physics
        td->fNphys++;
      }
      track.fPstep -= track.fSnext;
      track.fStep += track.fSnext;
      track.fSafety -= track.fSnext;
      if (track.fSafety < 0.)
        track.fSafety = 0;
      track.fXpos += track.fSnext * track.fXdir;
      track.fYpos += track.fSnext * track.fYdir;
      track.fZpos += track.fSnext * track.fZdir;
      // Update total number of steps
      td->fNsteps++;
      if (track.fSnext < 1.E-8) td->fNsmall++;
      track.fSnext = 0;
      MoveTrack(itr--, tracks, output);
    }
  }
  // Compact remaining tracks and move the removed oned to the output container
  // Check if tracking the remaining tracks can be postponed
  ntracks = tracks.size();
  action = PostponedAction(ntracks);
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
  double *steps = td->GetDblArray(ntracks);
  for (itr = 0; itr < fNtracks; itr++) {
    GeantTrack &track = *tracks[itr];
    lmax = SafeLength(track, eps);
    lmax = Math::Max<double>(lmax, track.fSafety);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = (track.fBoundary) ? 
                Math::Min<double>(lmax, Math::Max<double>(track.fSnext, 1.E-4)) 
              : Math::Min<double>(lmax, track.fPstep);
  }
  // Propagate the vector of tracks
  PropagateInVolume(tracks, ntracks, steps, td);
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
  for (int itr=0; itr<tracks.size(); ++itr) {
    track = *tracks[itr];
    if (track.fStatus == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      MoveTrack(itr--, tracks, output);
    }
  }
  // Select tracks which have not reached yet the boundary
  
  
  // Select tracks that are in flight or were propagated to boundary with
  // steps bigger than safety. Keep these tracks at the beginning of the
  // vector while moving the others to the end
  ntracks = GetNtracks();
  td->fNsteps += ntracks;

#if (defined(USE_VECGEOM_NAVIGATOR) && defined(VECTORIZED_GEOMETRY))
  icrossed += CheckSameLocation(tracks, ntracks, td);
#else
  for (itr = 0; itr < ntracks; itr++) {
    GeantTrack &track = *tracks[itr];
    int crossed = CheckSameLocationSingle(track, td);
    if (crossed) MoveTrack(itr--, tracks, output);
    icrossed += crossed;
  }
#endif
    
#ifdef BUG_HUNT
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateTracks");
#endif
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo::PropagateSingleTrack(TrackVec_t &tracks, int &itr, GeantTaskData *td, int stage) {
  // Propagate the track with its selected steps, starting from a given stage.

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
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep,
              "PropagateSingle", itr);
#endif

  TrackVec_t &output = *td->fTransported;
  GeantTrack &track = *tracks[itr];
  ComputeTransportLengthSingle(track, td);

#ifdef BUG_HUNT
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterCompTranspLenSingle");
#endif
  // Mark dead tracks for copy/removal
  if (track.fSnext < 0) {
    Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", track.fParticle);
    track.Print();
    track.fStatus = kKilled;
  }
  if (track.fStatus == kKilled) {
    MoveTrack(itr--, tracks, output); // next track at same index after move
    return 0;
  }
  // Stage 0: straight propagation
  if (stage == 0) {
    if (track.fCharge == 0 || bmag < 1.E-10) {
      // Do straight propagation to physics process or boundary
      if (track.fBoundary) {
        if (track.fNextpath->IsOutside())
          track.fStatus = kExitingSetup;
        else
          track.fStatus = kBoundary;
        icrossed++;
      } else {
        fStatusV[itr] = kPhysics;
        // Update number of steps to physics
        td->fNphys++;
      }
      track.fPstep -= track.fSnext;
      track.fStep += track.fSnext;
      track.fSafety -= track.fSnext;
      if (track.fSafety < 0.)
        track.fSafety = 0;
      track.fXpos += track.fSnext * track.fXdir;
      track.fYpos += track.fSnext * track.fYdir;
      track.fZpos += track.fSnext * track.fZdir;
      // Update total number of steps
      td->fNsteps++;
      if (track.fSnext < 1.E-8) td->fNsmall++;
      track.fSnext = 0;
      MoveTrack(itr--, tracks, output);

#ifdef BUG_HUNT
      BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateSingleNeutral",
                  track.fParticle);
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
    lmax = SafeLength(track, eps);
    lmax = Math::Max<double>(lmax, track.fSafety);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if frombdr=1) or the proposed  physics
    // step (frombdr=0)
    step = (track.fBoundary) ? 
             Math::Min<double>(lmax, Math::Max<double>(track.fSnext, 1.E-4)) 
           : Math::Min<double>(lmax, track.fPstep);
    // Propagate in magnetic field
    PropagateInVolumeSingle(track, step, td);
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
    if (track.fStatus == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      MoveTrack(itr--, tracks, output);
      return 0;
    }

    // Check if boundary has been crossed
    icrossed = CheckSameLocationSingle(track, td);
    if (icrossed) MoveTrack(itr--, tracks, output);
  }
#ifdef BUG_HUNT
  BreakOnStep(prop->fDebugEvt, prop->fDebugTrk, prop->fDebugStp, prop->fDebugRep, "AfterPropagateSingle", itr);
#endif
  return icrossed;
}

//______________________________________________________________________________
GEANT_CUDA_BOTH_CODE
int GeantTrackGeo::PropagateTracksScalar(TrackVec_t &tracks,
                                         GeantTaskData *td,
                                         int stage) {

  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  int icrossed = 0;
  for (int itr = 0; itr < tracks; ++itr) {
    icrossed += PropagateSingleTrack(tracks, itr, td, stage);
  }
  return icrossed;
}

//______________________________________________________________________________
int GeantTrackGeo::PostponeTracks(TrackVec_t &input, TrackVec_t &output) {
  // Postpone transport of remaining tracks and copy them to the output.
  int npostponed = input.size();
  assert(output.size() >= npostponed);
  // Move content
  std::move(input.begin(), input.end(), std::back_inserter(output));
  // Clear data
  input.clear();
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
