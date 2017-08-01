#include "FieldPropagationHandler.h"

#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "ConstFieldHelixStepper.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "ScalarNavInterfaceVG.h"
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
#else
#include "ScalarNavInterfaceTGeo.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::FieldPropagationHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::~FieldPropagationHandler()
{
// Destructor
}


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Scalar geometry length computation. The track is moved into the output basket.
  // Step selection
  double step, lmax;
  const double eps = 1.E-2; //
  const double bmag = fPropagator->fConfig->fBmag;

  // We use the track sagitta to estimate the "bending" error,
  // i.e. what is the propagated length for which the track deviation in
  // magnetic field with respect to straight propagation is less than epsilon.
  // Take the maximum between the safety and the "bending" safety
  lmax = SafeLength(*track, bmag, eps);
  lmax = vecCore::math::Max<double>(lmax, track->fSafety);
  // Select step to propagate as the minimum among the "safe" step and:
  // the straight distance to boundary (if frombdr=1) or the proposed  physics
  // step (frombdr=0)
  step = (track->fBoundary) ?
           vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track->fSnext, 1.E-4))
         : vecCore::math::Min<double>(lmax, track->fPstep);
  // Propagate in magnetic field
  PropagateInVolume(*track, step, td);
  //Update number of partial steps propagated in field
  td->fNmag++;
  // Update time of flight and number of interaction lengths
//  track->fTime += track->TimeStep(track->fStep);
//  track->fNintLen -= track->fStep/track->fIntLen;

  // Set continuous processes stage as follow-up for tracks that reached the
  // physics process
  if (track->fStatus == kPhysics) {
    // Update number of steps to physics and total number of steps
    td->fNphys++;
    td->fNsteps++;
#ifdef USE_REAL_PHYSICS
//    track->SetStage(kAlongStepActionStage);
    track->SetStage(kPostStepActionStage);
#else
    track->SetStage(kContinuousProcStage);
#endif
    output.AddTrack(track);
    return;
  }

  // Crossing tracks continue to continuous processes, the rest have to
  // query again the geometry
  if (!IsSameLocation(*track, td)) {
    td->fNcross++;
    td->fNsteps++;
  } else {
    track->SetStage(kGeometryStepStage);
  }
  output.AddTrack(track);
}


//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::DoIt(Basket &input, Basket& output, GeantTaskData *td)
{
// Vector geometry length computation. The tracks are moved into the output basket.
  TrackVec_t &tracks = input.Tracks();
  double lmax;
  const double eps = 1.E-2; //
  const double bmag = fPropagator->fConfig->fBmag;

  int ntracks = tracks.size();
  double *steps = td->GetDblArray(ntracks);
  for (int itr = 0; itr < ntracks; itr++) {
    // Can this loop be vectorized?
    GeantTrack &track = *tracks[itr];
    lmax = SafeLength(track, bmag, eps);
    lmax = vecCore::math::Max<double>(lmax, track.fSafety);
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = (track.fBoundary) ?
                vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track.fSnext, 1.E-4))
              : vecCore::math::Min<double>(lmax, track.fPstep);
  }
  // Propagate the vector of tracks
  PropagateInVolume(input.Tracks(), steps, td);

  //Update number of partial steps propagated in field
  td->fNmag += ntracks;

  // Update time of flight and number of interaction lengths.
  // Check also if it makes sense to call the vector interfaces
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  int nvect = 0;
#endif
  for (auto track : tracks) {
//    track->fTime += track->TimeStep(track->fStep);
//    track->fNintLen -= track->fStep/track->fIntLen;
    if (track->fStatus == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      output.AddTrack(track);
      continue;
    }
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
    if (track->fSafety < 1.E-10 || track->fSnext < 1.E-10)
      nvect++;
#else
    // Vector treatment was not requested, so proceed with scalar
    if (!IsSameLocation(*track, td)) {
      td->fNcross++;
      td->fNsteps++;
    } else {
      track->SetStage(kGeometryStepStage);
    }
    output.AddTrack(track);
    continue;
#endif
  }
  // If vectorized treatment was requested and the remaining population is
  // large enough, continue with vectorized treatment
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  constexpr int kMinVecSize = 8; // this should be retrieved from elsewhere
  if (nvect < kMinVecSize) {
    for (auto track : tracks) {
      if (track->fStatus == kPhysics) continue;
      if (!IsSameLocation(*track, td)) {
        td->fNcross++;
#ifdef USE_REAL_PHYSICS
//            track->SetStage(kAlongStepActionStage);
            track->SetStage(kPostStepActionStage);
#else
            track->SetStage(kContinuousProcStage);
#endif
      } else {
        track->SetStage(kGeometryStepStage);
      }
      output.AddTrack(track);
      continue;
    }
    return;
  }
  // This part deals with vectorized treatment
  // Copy data to SOA and dispatch for vector mode
  GeantTrackGeo_v &track_geo = *td.fGeoTrack;
  for (auto track : tracks) {
    if (track.fStatus != kPhysics &&
        (track.fSafety < 1.E-10 || track.fSnext < 1.E-10))
      track_geo.AddTrack(*track);
  }
  bool *same = td->GetBoolArray(nvect);
  NavigationState *tmpstate = td->GetPath();
  VectorNavInterface::NavIsSameLocation(nvect,
                   track_geo.fXposV, track_geo.fYposV, track_geo.fZposV,
                   track_geo.fXdirV, track_geo.fYdirV, track_geo.fZdirV,
                   (const VolumePath_t**)fPathV, fNextpathV, same, tmpstate);
  track_geo.UpdateOriginalTracks();
  for (itr = 0; itr < nsel; itr++) {
    GeantTrack *track = track_geo.fOriginalV[itr];
    if (!same[itr]) {
      td->fNcross++;
      td->fNsteps++;
      track->fBoundary = true;
      track->fStatus = kBoundary;
      if (track->fNextpath->IsOutside())
        track->fStatus = kExitingSetup;
      if (track->fStep < 1.E-8) td->fNsmall++;
    } else {
      track->fBoundary = false;
      track->SetStage(kGeometryStepStage);
    }
    output.AddTrack(track);
  }
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PropagateInVolume(TrackVec_t &tracks, const double *crtstep,
                                                GeantTaskData *td)
{
// THIS IS THE VECTORIZED IMPLEMENTATION PLACEHOLDER FOR MAGNETIC FIELD PROPAGATION.
// Now implemented just as a loop
  int ntracks = tracks.size();
  for (int itr=0; itr<ntracks; ++itr)
    PropagateInVolume(*tracks[itr], crtstep[itr], td);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PropagateInVolume(GeantTrack &track, double crtstep, GeantTaskData * td)
{
// Single track propagation in a volume. The method is to be called
// only with  charged tracks in magnetic field.The method decreases the fPstepV
// fSafetyV and fSnextV with the propagated values while increasing the fStepV.
// The status and boundary flags are set according to which gets hit first:
// - physics step (bdr=0)
// - safety step (bdr=0)
// - snext step (bdr=1)
   bool useRungeKutta;
   const double bmag = td->fPropagator->fConfig->fBmag;
   useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;

// #ifdef RUNGE_KUTTA
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
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
  if (track.fPstep < 1.E-10) {
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
  track.fStep += crtstep;
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
#ifndef VECCORE_CUDA
     fieldPropagator->DoStep(Position,    Direction,    track.fCharge, track.fP, crtstep,
                             PositionNew, DirectionNew);
#endif
  } else {
     // Old - constant field
     ConstBzFieldHelixStepper stepper(bmag);
     stepper.DoStep<ThreeVector,double,int>(Position,    Direction,    track.fCharge, track.fP, crtstep,
                                         PositionNew, DirectionNew);
  }

  track.fXpos = PositionNew.x();
  track.fYpos = PositionNew.y();
  track.fZpos = PositionNew.z();

  //  maybe normalize direction here  // vecCore::math::Normalize(dirnew);
  DirectionNew = DirectionNew.Unit();
  track.fXdir = DirectionNew.x();
  track.fYdir = DirectionNew.y();
  track.fZdir = DirectionNew.z();

#if 0
  ThreeVector SimplePosition = Position + crtstep * Direction;
  // double diffpos2 = (PositionNew - Position).Mag2();
  double diffpos2 = (PositionNew - SimplePosition).Mag2();
  //   -- if (vecCore::math::Sqrt(diffpos)>0.01*crtstep) {
  const double drift= 0.01*crtstep;
  if ( diffpos2>drift*drift ){
      double diffpos= vecCore::math::Sqrt(diffpos2);
      // Geant::Print("PropagateInVolumeSingle","relative difference in pos = %g", diffpos/crtstep);
      Geant::Print("PropagateInVolumeSingle","difference in pos = %g (abs) %g (relative) , step= %g",
                   diffpos, diffpos/crtstep, crtstep);
  }
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool FieldPropagationHandler::IsSameLocation(GeantTrack &track, GeantTaskData *td) {
// Query geometry if the location has changed for a track
// Returns number of tracks crossing the boundary (0 or 1)

  if (track.fSafety > 1.E-10 && track.fSnext > 1.E-10) {
    // Track stays in the same volume
    track.fBoundary = false;
    return true;
  }

  // Track may have crossed, check it
  bool same;
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::NavigationState *tmpstate = td->GetPath();
  ScalarNavInterfaceVGM::NavIsSameLocation(track, same, tmpstate);
#else
// ROOT navigation
  ScalarNavInterfaceTGeo::NavIsSameLocation(track, same);
#endif // USE_VECGEOM_NAVIGATOR
  if (same) {
    track.fBoundary = false;
    return true;
  }

  track.fBoundary = true;
  track.fStatus = kBoundary;
  if (track.NextPath()->IsOutside())
    track.fStatus = kExitingSetup;
  if (track.fStep < 1.E-8) td->fNsmall++;
  return false;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
