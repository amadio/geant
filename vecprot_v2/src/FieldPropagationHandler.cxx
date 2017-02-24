#include "FieldPropagationHandler.h"

#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "ConstFieldHelixStepper.h"

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
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  const double bmag = gPropagator_fConfig->fBmag;
#else
  const double bmag = fPropagator->fConfig->fBmag;
#endif
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
  // Update time of flight and number of interaction lengths
  track->fTime += track->TimeStep(track->fStep);
  track->fNintLen -= track->fStep/track->fIntLen;
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
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  const double bmag = gPropagator_fConfig->fBmag;
#else
  const double bmag = fPropagator->fConfig->fBmag;
#endif
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

  // Update time of flight and number of interaction lengths
  for (auto track : tracks) {
    track->fTime += track->TimeStep(track->fStep);
    track->fNintLen -= track->fStep/track->fIntLen;  
  }
  
  // Copy tracks to output
#ifndef VECCORE_CUDA
  std::move(tracks.begin(), tracks.end(), std::back_inserter(output.Tracks()));
#else
  for (auto track : tracks) output->AddTrack(track);
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
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
   const double bmag = gPropagator_fConfig->fBmag;
   constexpr auto gPropagator_fUseRK = false; // Temporary work-around until actual implementation ..
   useRungeKutta= gPropagator_fUseRK;   //  Something like this is needed - TBD
#else
   const double bmag = td->fPropagator->fConfig->fBmag;
   useRungeKutta= td->fPropagator->fConfig->fUseRungeKutta;
#endif

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

} // GEANT_IMPL_NAMESPACE
} // Geant
