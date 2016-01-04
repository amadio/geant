#include "TransportManager.h"

#include "globals.h"
#include "Geant/Error.h"
#include <execinfo.h>

#include "GeantTrackGeo.h"
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

#include "Basket.h"
#include "GeantTaskData.h"
#include "ConstBzFieldHelixStepper.h"
#include "ConstVecFieldHelixStepper.h"
#include "GeantScheduler.h"

// #ifdef  RUNGE_KUTTA
#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"

#include "GUVField.h"
#include "FieldLookup.h"

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
VECCORE_ATT_HOST_DEVICE
int TransportManager::CheckSameLocation(TrackVec_t &tracks,
                                         int ntracks,
                                         GeantTaskData *td) {
// Query geometry if the location has changed for a vector of particles
// Returns number of tracks crossing the boundary.

  TrackVec_t &output = td->fTransported1;
  int icrossed = 0;
  int nsel = ntracks;
  for (int itr = 0; itr < ntracks; itr++) {
    GeantTrack &track = *tracks[itr];
    if (track.GetSafety() > 1.E-10 && track.GetSnext() > 1.E-10) {
      // Track stays in the same volume
      track.SetBoundary(false);
      nsel--;
    }
  }
  if (!nsel) return 0;
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  constexpr int kMinVecSize = 4; // this should be retrieved from elsewhere
  if (nsel < kMinVecSize) {
    // No efficient vectorization possible, do scalar treatment
#endif
    // The scalar treatment is always done in case of non-vector compilation
    for (unsigned int itr = 0; itr < tracks.size(); itr++) {
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
  GeantTrackGeo_v &track_geo = *td.fGeoTrack;
  track_geo.Clear();
  for (int itr = 0; itr < tracks.size(); itr++) {
    GeantTrack &track = *tracks[itr];
    if (track.GetSafety() < 1.E-10 || track.GetSnext() < 1.E-10)
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
    if (track.Boundary()) MoveTrack(itr--, tracks, output);
  }

  return icrossed;
#endif   // vector mode
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TransportManager::CheckSameLocationSingle(GeantTrack &track,
                                         GeantTaskData *td) {
// Query geometry if the location has changed for a track
// Returns number of tracks crossing the boundary (0 or 1)

  if (track.GetSafety() > 1.E-10 && track.GetSnext() > 1.E-10) {
    // Track stays in the same volume
    track.SetBoundary(false);
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
    track.SetBoundary(false);
    return 0;
  }

  track.SetBoundary(true);
  track.SetStatus(kBoundary);
  if (track.NextPath()->IsOutside())
    track.SetStatus(kExitingSetup);
  if (track.GetStep() < 1.E-8) td->fNsmall++;
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TransportManager::ComputeTransportLength(TrackVec_t &tracks,
                                              int ntracks,
                                              GeantTaskData *td) {
// Vector version for proposing the geometry step. All tracks have to be
// in the same volume. This may still fall back on the scalar implementation
// in case the vector size is too small.
#ifdef USE_VECGEOM_NAVIGATOR
//#define VECTORIZED_GEOMETRY
#ifdef VECTORIZED_GEOMETRY
  constexpr int kMinVecSize = 4; // this should be retrieved from elsewhere
  // We ignore tracks for which the current safety allows for the proposed step
  // We need to count if the remaining tracks may form a vector
  int nsel = 0;
  for (auto track : tracks)
    if (track->GetSafety() < track->GetPstep()) nsel++;
  if (nsel < kMinVecSize) {
    // Just process all the tracks in scalar mode
    for (auto track : tracks)
      ComputeTransportLengthSingle(*track, td);
  } else {
    // Copy interesting tracks to SOA and process vectorized.
    // This is overhead but cannot be avoided.
    GeantTrackGeo_v &track_geo = *td->fGeoTrack;
    track_geo.Clear();
    int i = 0;
    // This selection should happen in the propagation stage
    for (auto track : tracks) {
      if (track->GetSafety() < track->GetPstep()) {
        track_geo.AddTrack(*track);
        td->fPathV[i] = track->Path();
        td->fNextpathV[i] = track->NextPath();
        i++;
      }
    }
    // The vectorized SOA call
    VectorNavInterface::NavFindNextBoundaryAndStep(nsel, track_geo.fPstepV,
                 track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
                 track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
                 (const VolumePath_t **)td->fPathV, td->fNextpathV,
                 track_geo.fSnextV, track_geo.fSafetyV, track_geo.fBoundaryV);

    // Update original tracks
    track_geo.UpdateOriginalTracks();
    // Update number of calls to geometry (1 vector call + ntail scalar calls)
    td->fNsnext += nsel;
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
    if ((track->NextPath()->IsOutside() && track->GetSnext() < 1.E-6) || track->GetSnext() > 1.E19)
      track->SetStatus(kExitingSetup);
  }
#else
  // TGeo implementation fall back on looped version
  (void)ntracks;
  for (auto track : tracks) {
    ComputeTransportLengthSingle(*track, td);
  }
#endif // USE_VECGEOM_NAVIGATOR
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TransportManager::ComputeTransportLengthSingle(GeantTrack &track, GeantTaskData *td) {
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
  if (track.NextPath()->IsOutside() || track.GetSnext() > 1.E19)
    track.SetStatus(kExitingSetup);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TransportManager::PropagateInVolume(TrackVec_t &tracks,
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
  for (int itr=0; itr<ntracks; ++itr) {
    GeantTrack &track = *tracks[itr];
    PropagateInVolumeSingle(track, crtstep[itr++], td);
  }
}


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TransportManager::PropagateInVolumeSingle(GeantTrack &track, double crtstep, GeantTaskData * td) {
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
   // using Geant::ConstBzFieldHelixStepper;
   // using Geant::cxx::ConstBzFieldHelixStepper;   
   
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
   constexpr auto gPropagator_fUseRK = false; // Temporary work-around until actual implementation ..
   useRungeKutta= gPropagator_fUseRK;   //  Something like this is needed - TBD
#else
   const double bmag = td->fPropagator->fConfig->fBmag;
   useRungeKutta= td->fPropagator->fConfig->fUseRungeKutta;
#endif

   // static unsigned long icount= 0;
   // if( icount++ < 2 )  std::cout << " PropagateInVolumeSingle: useRungeKutta= " << useRungeKutta << std::endl;

// #ifdef RUNGE_KUTTA
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
   GUFieldPropagator *fieldPropagator = nullptr;
   if( useRungeKutta ){
      // Initialize for the current thread -- move to GeantPropagator::Initialize() or per thread Init method
      static GUFieldPropagatorPool* fieldPropPool= GUFieldPropagatorPool::Instance();
      assert( fieldPropPool );

      fieldPropagator = fieldPropPool->GetPropagator(td->fTid);
      assert( fieldPropagator );
      td->fFieldPropagator= fieldPropagator;
   }
#endif

  // Reset relevant variables
  track.SetStatus(kInFlight);
  track.DecreasePstep(crtstep);
  if (track.GetPstep() < 1.E-10) {
    track.SetPstep(0);
    track.SetStatus(kPhysics);
  }
  track.DecreaseSafety(crtstep);
  if (track.GetSafety() < 1.E-10)
    track.SetSafety(0);
  track.DecreaseSnext(crtstep);
  if (track.GetSnext() < 1.E-10) {
    track.SetSnext(0);
    if (track.Boundary()) {
      track.SetStatus(kBoundary);
    }
  }
  track.IncreaseStep(crtstep);
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
// alternative code with lean stepper would be:
// ( stepper header has to be included )

  using ThreeVector = vecgeom::Vector3D<double>;
  ThreeVector Position(track.X(), track.Y(), track.Z());
  ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

  double curvaturePlus= fabs(GeantTrack::kB2C * track.fCharge * bmag) / (track.fP + 1.0e-30);  // norm for step
  // 'Curvature' along the full track - not just in the plane perpendicular to the B-field vector

  constexpr double numRadiansMax= 10.0;   //  Too large an angle - many RK steps.  Potential change -> 2.0*PI;
  constexpr double numRadiansMin= 0.05;   //  Very small an angle - helix is adequate.  TBC: Use average B-field value?
      //  A track turning more than 10 radians will be treated approximately
  const double angle= crtstep * curvaturePlus;
  bool mediumAngle = ( numRadiansMin < angle ) && ( angle < numRadiansMax );
  useRungeKutta = useRungeKutta && mediumAngle;
  
  // double BfieldInitial[3], bmag= 0.0;
  // FieldLookup::GetFieldValue(td, Position, BfieldInitial, &bmag);
  // printf("TransportMgr::PropagateInVolumeSingle> Curvature= %8.4g  CurvPlus= %8.4g  step= %f  Bmag=%8.4g  momentum mag=%f  angle= %g\n",
  //       Curvature(td, i), curvaturePlus, crtstep, bmag, track.fP, angle );

  // Option: use RK as fall back - until 'General Helix' is robust
  // bool dominantBz =  std::fabs( std::fabs(BfieldInitial[1]) )
  //       > 1e6 *
  //   std::max( std::fabs( BfieldInitial[0]), std::fabs(BfieldInitial[1]) );
  // if( !dominantBz )
  //   useRungeKutta = true;
  // int propagationType= 0;
  
  if( useRungeKutta ) {
#ifndef VECCORE_CUDA
     fieldPropagator->DoStep(Position,    Direction,    track.Charge(), track.P(), crtstep,
                             PositionNew, DirectionNew);
     // crtstep = 1.0e-4;   printf( "Setting crtstep = %f -- for DEBUGing ONLY.", crtstep );
     // propagationType= 1;

     // CheckDirection(DirectionNew);
     
     /**
     const bool fCheckingStep= false;
     if( fCheckingStep ) {
        const double epsDiff= 2.0e-3; // bool verbDiff= true );
        StepChecker EndChecker( epsDiff, epsDiff * crtstep, true );
        vecgeom::Vector3D<double> Bfield( BfieldInitial[0], BfieldInitial[1], BfieldInitial[2] );
        EndChecker.CheckStep( Position, Direction, track.fCharge, track.fP, crtstep,
                              PositionNewRK, DirectionNewRK, Bfield );
     }  
     **/   
#endif
  } else {
     double BfieldInitial[3], bmag= 0.0;
     FieldLookup::GetFieldValue(td, Position, BfieldInitial, &bmag);
     double Bx= BfieldInitial[0], By= BfieldInitial[1], Bz= BfieldInitial[2];
     if ( std::fabs( Bz ) > 1.0e6 * std::max( std::fabs(Bx), std::fabs(By) ) )
     {
        ConstBzFieldHelixStepper stepper( Bz );
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), crtstep,
                                               PositionNew, DirectionNew);
        // propagationType= 2;
     } else {
        ConstVecFieldHelixStepper stepper( BfieldInitial );
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), crtstep,
                                         PositionNew, DirectionNew);
        // propagationType= 3;        
     }     
  }

  track.SetPosition(PositionNew);

  //  maybe normalize direction here  // Math::Normalize(dirnew);
  DirectionNew = DirectionNew.Unit();
  track.SetDirection(DirectionNew);

#if 0
  ThreeVector SimplePosition = Position + crtstep * Direction;
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
int TransportManager::PropagateTracks(TrackVec_t &tracks, GeantTaskData *td) {
  // Propagate the ntracks in the current volume with their physics steps (already
  // computed)
  // Vectors are pushed downstream when efficient.
  TrackVec_t &output = td->fTransported1;
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
  GeantPropagator *prop = td->fPropagator;  
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "PropagateTracks");
#endif
  ComputeTransportLength(tracks, ntracks, td);
//         Printf("====== After ComputeTransportLength:");
//         PrintTracks();
#ifdef BUG_HUNT
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterCompTransLen");
#endif

  int itr = 0;
  int icrossed = 0;
  double lmax;
  const double eps = 1.E-2; // 100 micron
  const double bmag = td->fBfieldMag;

  // Remove dead tracks, propagate neutrals
  for (unsigned int itr=0; itr<tracks.size(); ++itr) {
    GeantTrack &track = *tracks[itr];
    // Move dead to output
    if (track.GetSnext() < 0) {
      Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", track.Particle());
      track.Print("");
      track.SetStatus(kKilled);
    }
    if (track.Status() == kKilled) {
      MoveTrack(itr--, tracks, output); // next track at same index after move
      continue;
    }
    // Propagate straight tracks to the precomputed location and update state,
    // ten move them to output
    if (track.Charge() == 0 || bmag < 1.E-10) {
      // Do straight propagation to physics process or boundary
      if (track.Boundary()) {
        if (track.NextPath()->IsOutside())
          track.SetStatus(kExitingSetup);
        else
          track.SetStatus(kBoundary);
        icrossed++;
      } else {
        track.SetStatus(kPhysics);
        // Update number of steps to physics
        td->fNphys++;
      }
      track.MakeStep(track.GetSnext());
      if (track.GetSafety() < 0.)
        track.SetSafety(0);
      // Update total number of steps
      td->fNsteps++;
      if (track.GetSnext() < 1.E-8) td->fNsmall++;
      track.SetSnext(0);
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
    icrossed += TransportManager::PropagateTracksScalar(tracks, td, 1);
    return icrossed;
  case kPostpone:
    PostponeTracks(tracks, output);
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
  double *steps = td->GetDblArray(ntracks);
  for (itr = 0; itr < ntracks; itr++) {
    GeantTrack &track = *tracks[itr];
    lmax = SafeLength(track, bmag, eps);
    lmax = Math::Max<double>(lmax, track.GetSafety());
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = (track.Boundary()) ?
                Math::Min<double>(lmax, Math::Max<double>(track.GetSnext(), 1.E-4))
              : Math::Min<double>(lmax, track.GetPstep());
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
  for (unsigned int itr=0; itr<tracks.size(); ++itr) {
    GeantTrack &track = *tracks[itr];
    if (track.Status() == kPhysics) {
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
  ntracks = tracks.size();
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
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateTracks");
#endif
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TransportManager::PropagateSingleTrack(GeantTrack *track, Basket *output, GeantTaskData *td, int stage)
{
  // Propagate the track with its selected steps, starting from a given stage.
  GeantPropagator *prop = td->fPropagator;
  int icrossed = 0;
  double step, lmax;
  const double eps = 1.E-2; // 1 micron
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  const double bmag = gPropagator_fConfig->fBmag;
#else
  const double bmag = prop->fConfig->fBmag;
#endif
// Compute transport length in geometry, limited by the physics step
  ComputeTransportLengthSingle(*track, td);
  // Mark dead tracks for copy/removal
  if (track->GetSnext() < 0) {
    Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", track->Particle());
    track->Print("");
    track->SetStatus(kKilled);
  }
  if (track->Status() == kKilled) {
    output->AddTrack(track);
    return 0;
  }
  // Stage 0: straight propagation for neutrals
  if (stage == 0) {
    if (track->Charge() == 0 || bmag < 1.E-10) {
      // Do straight propagation to physics process or boundary
      if (track->Boundary()) {
        if (track->NextPath()->IsOutside())
          track->SetStatus(kExitingSetup);
        else
          track->SetStatus(kBoundary);
        icrossed++;
      } else {
        track->SetStatus(kPhysics);
        // Update number of steps to physics
        td->fNphys++;
      }
      track->MakeStep(track->GetSnext());
      if (track->GetSafety() < 0.)
        track->SetSafety(0);
      // Update total number of steps
      td->fNsteps++;
      if (track->GetSnext() < 1.E-8) td->fNsmall++;
      track->SetSnext(0);
      output->AddTrack(track);
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
    lmax = SafeLength(*track, bmag, eps);
    lmax = Math::Max<double>(lmax, track->GetSafety());
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if frombdr=1) or the proposed  physics
    // step (frombdr=0)
    step = (track->Boundary()) ?
             Math::Min<double>(lmax, Math::Max<double>(track->GetSnext(), 1.E-4))
           : Math::Min<double>(lmax, track->GetPstep());
    // Propagate in magnetic field
    PropagateInVolumeSingle(*track, step, td);
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
    if (track->Status() == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      output->AddTrack(track);
      return 0;
    }

    // Check if boundary has been crossed
    icrossed = CheckSameLocationSingle(*track, td);
    if (icrossed) output->AddTrack(track);
  }
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TransportManager::PropagateSingleTrack(TrackVec_t &tracks, int &itr, GeantTaskData *td, int stage) {
  // Propagate the track with its selected steps, starting from a given stage.

  int icrossed = 0;
  double step, lmax;
  const double eps = 1.E-2; // 1 micron

  double Bfield[3], bmag= 0.0;
  // const double bmag = td->fBfieldMag;
  
// Compute transport length in geometry, limited by the physics step

#ifdef BUG_HUNT
  GeantPropagator *prop = td->fPropagator;
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep,
              "PropagateSingle", itr);
#endif

  TrackVec_t &output = td->fTransported1;
  GeantTrack &track = *tracks[itr];
  ComputeTransportLengthSingle(track, td);

#ifdef BUG_HUNT
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterCompTranspLenSingle");
#endif
  // Mark dead tracks for copy/removal
  if (track.GetSnext() < 0) {
    Error("ComputeTransportLength", "Track %d cannot cross boundary and has to be killed", track.Particle());
    track.Print("");
    track.SetStatus(kKilled);
  }
  if (track.Status() == kKilled) {
    MoveTrack(itr--, tracks, output); // next track at same index after move
    return 0;
  }
  // Stage 0: straight propagation
  if (stage == 0) {
    bool neutral = (track.fCharge == 0);
    if( !neutral ) {
       // printf( " PropagateSingleTrack> getting Field. Charge= %3d ", track.fCharge );
       vecgeom::Vector3D<double> Position( tracks[itr]->X(), tracks[itr]->Y(), tracks[itr]->Z() );
       FieldLookup::GetFieldValue(td, Position, Bfield, &bmag);       
       // if( bmag < 1.E-10) { printf("TransportMgr::TrSnglTrk> Tiny field - mag = %f\n", bmag); }
    }
    if ( neutral || bmag < 1.E-10) {       
      // Do straight propagation to physics process or boundary
      if (track.Boundary()) {
        if (track.NextPath()->IsOutside())
          track.SetStatus(kExitingSetup);
        else
          track.SetStatus(kBoundary);
        icrossed++;
      } else {
        track.SetStatus(kPhysics);
        // Update number of steps to physics
        td->fNphys++;
      }
      track.MakeStep(track.GetSnext());
      if (track.GetSafety() < 0.)
        track.SetSafety(0);
      // Update total number of steps
      td->fNsteps++;
      if (track.GetSnext() < 1.E-8) td->fNsmall++;
      track.SetSnext(0);
      MoveTrack(itr--, tracks, output);

#ifdef BUG_HUNT
      BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateSingleNeutral",
                  track.Particle());
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
    lmax = SafeLength(track, bmag, eps);
    lmax = Math::Max<double>(lmax, track.GetSafety());
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if frombdr=1) or the proposed  physics
    // step (frombdr=0)
    step = (track.Boundary()) ?
             Math::Min<double>(lmax, Math::Max<double>(track.GetSnext(), 1.E-4))
           : Math::Min<double>(lmax, track.GetPstep());
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
    if (track.Status() == kPhysics) {
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
  BreakOnStep(tracks, prop->fConfig->fDebugEvt, prop->fConfig->fDebugTrk, prop->fConfig->fDebugStp, prop->fConfig->fDebugRep, "AfterPropagateSingle", itr);
#endif
  return icrossed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int TransportManager::PropagateTracksScalar(TrackVec_t &tracks,
                                         GeantTaskData *td,
                                         int stage) {

  // Propagate the tracks with their selected steps in a single loop,
  // starting from a given stage.

  int icrossed = 0;
  for (int itr = 0; itr < (int)tracks.size(); ++itr) {
    icrossed += PropagateSingleTrack(tracks, itr, td, stage);
  }
  return icrossed;
}

//______________________________________________________________________________
int TransportManager::PostponeTracks(TrackVec_t &input, TrackVec_t &output) {
  // Postpone transport of remaining tracks and copy them to the output.
  auto npostponed = input.size();
  assert(output.size() >= npostponed);
  // Move content
  std::move(input.begin(), input.end(), std::back_inserter(output));
  // Clear data
  input.clear();
  return npostponed;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool TransportManager::BreakOnStep(TrackVec_t &tracks, int evt, int trk, int stp, int nsteps, const char *msg, int itr) {
  // Return true if container has a track with a given number doing a given step from a given event
  // Debugging purpose
  int start = 0;
  int end = tracks.size();
  bool has_it = false;
  if (itr >= 0) {
    start = itr;
    end = itr + 1;
  }
  for (itr = start; itr < end; ++itr) {
    GeantTrack &track = *tracks[itr];
    if ((track.Particle() == trk) && (track.Event() == evt) &&
        ((track.GetNsteps() >= stp) && (track.GetNsteps() < stp + nsteps))) {
      has_it = true;
#ifndef VECCORE_CUDA
      track.Print(msg);
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

} // Geant
