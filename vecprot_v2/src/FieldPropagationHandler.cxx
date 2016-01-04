#include "FieldPropagationHandler.h"

#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "ConstBzFieldHelixStepper.h"
#include "ConstVecFieldHelixStepper.h"

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
  lmax = vecCore::math::Max<double>(lmax, track->GetSafety());
  // Select step to propagate as the minimum among the "safe" step and:
  // the straight distance to boundary (if frombdr=1) or the proposed  physics
  // step (frombdr=0)
  step = (track->Boundary()) ?
           vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track->GetSnext(), 1.E-4))
         : vecCore::math::Min<double>(lmax, track->GetPstep());
  // Propagate in magnetic field
  PropagateInVolume(*track, step, td);
  //Update number of partial steps propagated in field
  td->fNmag++;

  // Set continuous processes stage as follow-up for tracks that reached the
  // physics process
  if (track->Status() == kPhysics) {
    // Update number of steps to physics and total number of steps
    td->fNphys++;
    td->fNsteps++;
#ifdef USE_REAL_PHYSICS
//    track->SetStage(kAlongStepActionStage);
    track->SetStage(kPostPropagationStage);
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
    lmax = vecCore::math::Max<double>(lmax, track.GetSafety());
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = (track.Boundary()) ?
                vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track.GetSnext(), 1.E-4))
              : vecCore::math::Min<double>(lmax, track.GetPstep());
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
    if (track->Status() == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;
      td->fNsteps++;
      output.AddTrack(track);
      continue;
    }
#if (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
    if (track->GetSafety() < 1.E-10 || track->GetSnext() < 1.E-10)
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
      if (track->Status() == kPhysics) continue;
      if (!IsSameLocation(*track, td)) {
        td->fNcross++;
#ifdef USE_REAL_PHYSICS
            track->SetStage(kPostPropagationStage);
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
    if (track.Status() != kPhysics &&
        (track.GetSafety() < 1.E-10 || track.GetSnext() < 1.E-10))
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
      track->SetBoundary(true);
      track->SetStatus(kBoundary);
      if (track->NextPath()->IsOutside())
        track->SetStatus(kExitingSetup);
      if (track->GetStep() < 1.E-8) td->fNsmall++;
    } else {
      track->SetBoundary(false);
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
   using ThreeVector = vecgeom::Vector3D<double>;  
   bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;
   // double bmag = td->fPropagator->fConfig->fBmag;
   double BfieldInitial[3], bmag;
   ThreeVector Position(track.X(), track.Y(), track.Z());
   FieldLookup::GetFieldValue(td, Position, BfieldInitial, &bmag);   

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
   GUFieldPropagator *fieldPropagator = nullptr;
   if( useRungeKutta ){
      fieldPropagator = td->fFieldPropagator;
      assert( fieldPropagator );
   }
   // GUFieldPropagator *fieldPropagator = useRungeKutta ? td->fFieldPropagator : nullptr;   
#endif

  // Reset relevant variables
  track.SetStatus(kInFlight);
  double pstep = track.GetPstep() - crtstep;
  if (pstep < 1.E-10) {
    pstep = 0;
    track.SetStatus(kPhysics);
  }
  track.SetPstep(pstep);
  double safety = track.GetSafety() - crtstep;
  if (safety < 1.E-10) safety = 0;
  track.SetSafety(safety);
  double snext = track.GetSnext() - crtstep;
  if (snext < 1.E-10) {
    snext = 0;
    if (track.Boundary())
      track.SetStatus(kBoundary);
  }
  track.SetSnext(snext);
  track.IncreaseStep(crtstep);
#ifdef USE_VECGEOM_NAVIGATOR
//  CheckLocationPathConsistency(i);
#endif
  double curvaturePlus= fabs(GeantTrack::kB2C * track.fCharge * bmag) / (track.fP + 1.0e-30);  // norm for step
  
  constexpr double numRadiansMax= 10.0; // Too large an angle - many RK steps.  Potential change -> 2.0*PI;
  constexpr double numRadiansMin= 0.05; // Very small an angle - helix is adequate.  TBC: Use average B-field value?
      //  A track turning more than 10 radians will be treated approximately
  const double angle= crtstep * curvaturePlus;
  bool mediumAngle = ( numRadiansMin < angle ) && ( angle < numRadiansMax );
  useRungeKutta = useRungeKutta && (mediumAngle);

  bool dominantBz =  std::fabs( std::fabs(BfieldInitial[2]) )
     > 1.e3 * std::max( std::fabs( BfieldInitial[0]), std::fabs(BfieldInitial[1]) );

#ifdef DEBUG_FIELD
  printf("--PropagateInVolume(Single): \n");
  printf("Curvature= %8.4g   CurvPlus= %8.4g  step= %f   Bmag=%8.4g   momentum mag=%f  angle= %g\n"
         Curvature(td, i), curvaturePlus, crtstep, bmag, fPV[i], angle );
#endif

  ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

#ifndef VECCORE_CUDA
  if( useRungeKutta ) {
     fieldPropagator->DoStep(Position,    Direction,    track.Charge(), track.P(), crtstep,
                             PositionNew, DirectionNew);
  }
  else
#endif
  {
     constexpr double toKiloGauss= 1.0e+14; // Converts to kilogauss -- i.e. 1 / Unit::kilogauss
                                            // Must agree with values in magneticfield/inc/Units.h
     double Bz = BfieldInitial[2] * toKiloGauss;
     if ( dominantBz ) {
        // Constant field in Z-direction
        ConstBzFieldHelixStepper stepper( Bz ); //
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,    track.Charge(), track.P(), crtstep,
                                               PositionNew, DirectionNew);
     } else {
        // Geant::
        ConstVecFieldHelixStepper stepper( BfieldInitial );
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), crtstep,                                           
                                               PositionNew, DirectionNew);        
     }
  }

  track.SetPosition(PositionNew);

  //  maybe normalize direction here  // vecCore::math::Normalize(dirnew);
  DirectionNew = DirectionNew.Unit();
  track.SetDirection(DirectionNew);

#ifdef REPORT_AND_CHECK
  double origMag= Direction.Mag();
  double oldMag= DirectionNew.Mag();
  double newMag= DirectionUnit.Mag();  
  Printf(" -- State after propagation in field:  Position= %f, %f, %f   Direction= %f, %f, %f  - mag original, integrated, integr-1.0, normed-1 = %10.8f %10.8f %7.2g %10.8f %7.2g",
         track.X(),  track.Y(),  track.Z(),
         track.Dx(),  track.Dy(),  track.Dz(),         
         origMag, oldMag, oldMag-1.0, newMag, newMag-1.0 );
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

  if (track.GetSafety() > 1.E-10 && track.GetSnext() > 1.E-10) {
    // Track stays in the same volume
    track.SetBoundary(false);
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
    track.SetBoundary(false);
    return true;
  }

  track.SetBoundary(true);
  track.SetStatus(kBoundary);
  if (track.NextPath()->IsOutside())
    track.SetStatus(kExitingSetup);
  if (track.GetStep() < 1.E-8) td->fNsmall++;
  return false;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
