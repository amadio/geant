#include "FieldPropagationHandler.h"

#include "FieldConfig.h"
#include "FieldLookup.h"

#include "GUFieldPropagatorPool.h"
#include "GUFieldPropagator.h"
#include "ConstBzFieldHelixStepper.h"
#include "ConstFieldHelixStepper.h"
#include "FieldTrack.h"

#include "GeantTrack.h"

#include "base/SOA3D.h"
// #include "SOA6D.h"
#include "Geant/VectorTypes.h"   // Defines Geant::Double_v etc

#include "WorkspaceForFieldPropagation.h"
#include "FlexIntegrationDriver.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
#include "ScalarNavInterfaceVG.h"
#include "ScalarNavInterfaceVGM.h"
#include "VectorNavInterface.h"
#else
#include "ScalarNavInterfaceTGeo.h"
#endif

using Double_v = Geant::Double_v;

// #define CHECK_VS_SCALAR  1

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

const double FieldPropagationHandler::gEpsDeflection = 1.E-2 * geant::cm; //Units
          
#ifdef USE_REAL_PHYSICS  
  auto stageAfterCrossing= kPostPropagationStage;
#else
  auto stageAfterCrossing= kContinuousProcStage;
#endif
          
#ifdef STATS_METHODS
static std::atomic<unsigned long> numRK      ,
                                  numHelixZ  ,
                                  numHelixGen,
                                  numTot;
#endif

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::FieldPropagationHandler(int threshold, GeantPropagator *propagator)
               : Handler(threshold, propagator)
{
// Default constructor
   // std::cout << " FieldPropagationHandler c-tor called:  threshold= " << threshold << std::endl;

#ifdef STATS_METHODS
   numTot      = 0;
   numRK       = 0;
   numHelixZ   = 0;
   numHelixGen = 0;
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::~FieldPropagationHandler()
{
// Destructor
}

//______________________________________________________________________________
GUFieldPropagator *
FieldPropagationHandler::Initialize(GeantTaskData * td)
{
  GUFieldPropagator *fieldPropagator = nullptr;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;

  if( useRungeKutta ){
     // Initialize for the current thread -- move to GeantPropagator::Initialize() or per thread Init method
     static GUFieldPropagatorPool* fieldPropPool= GUFieldPropagatorPool::Instance();
     assert( fieldPropPool );

     fieldPropagator = fieldPropPool->GetPropagator(td->fTid);
     assert( fieldPropagator );
     td->fFieldPropagator= fieldPropagator;
  }

  size_t basketSize = td->fPropagator->fConfig->fNperBasket;
  td->fSpace4FieldProp = new WorkspaceForFieldPropagation(basketSize);
#endif
  // For the moment return the field Propagator.
  // Once this method is called by the framework, this can be obtained from td->fFieldPropagator
  return fieldPropagator;
}

//______________________________________________________________________________
void FieldPropagationHandler::Cleanup(GeantTaskData * td)
{
   delete td->fSpace4FieldProp;
   td->fSpace4FieldProp = nullptr;
}

//______________________________________________________________________________
// Curvature for general field
VECCORE_ATT_HOST_DEVICE
double FieldPropagationHandler::Curvature(const GeantTrack  & track) const
{
  using ThreeVector_d = vecgeom::Vector3D<double>;
  constexpr double tiny = 1.E-30;
  constexpr double inv_kilogauss = 1.0 / geant::kilogauss;
  ThreeVector_d MagFld;
  double bmag= 0.0;

  ThreeVector_d Position(track.X(), track.Y(), track.Z());
  FieldLookup::GetFieldValue(Position, MagFld, bmag); // , td);
  bmag *= inv_kilogauss;

  //  Calculate transverse momentum 'Pt' for field 'B'
  //
  ThreeVector_d Momentum( track.Dx(), track.Dy(), track.Dz() );
  Momentum *= track.P();
  ThreeVector_d PtransB;  //  Transverse wrt direction of B
  double ratioOverFld = 0.0;
  if( bmag > 0 ) ratioOverFld = Momentum.Dot( MagFld ) / (bmag*bmag);
  PtransB = Momentum - ratioOverFld * MagFld ;
  double Pt_mag = PtransB.Mag();

  return fabs(GeantTrack::kB2C * track.Charge() * bmag / (Pt_mag + tiny));
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::DoIt(GeantTrack *track, Basket& output, GeantTaskData *td)
{
// Scalar geometry length computation. The track is moved into the output basket.
  // Step selection
  double step, lmax;
  const double epsDeflection = 1.E-2 * geant::cm;  // Units!
  
  // std::cout <<" FieldPropagationHandler::DoIt(*track) called for 1 ptrTrack." << std::endl;

  // We use the track sagitta to estimate the "bending" error,
  // i.e. what is the propagated length for which the track deviation in
  // magnetic field with respect to straight propagation is less than epsilon.
  // Take the maximum between the safety and the "bending" safety
  lmax = SafeLength( *track, /*gEps*/ epsDeflection);  
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
    track->SetStage(stageAfterCrossing); // Future: (kPostPropagationStage);    
  } else {
    // Crossing tracks continue to continuous processes, the rest have to
    // query again the geometry
    if (!IsSameLocation(*track, td)) {
      td->fNcross++;
      td->fNsteps++;
    } else {
      track->SetStage(kGeometryStepStage);
    }
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
  // const double gEpsDeflection = 1.E-2 * geant::cm;  // Units!
  
  // using minD= vecCore::math::Min<double>;
  // using maxD= vecCore::math::Max<double>;

  int ntracks = tracks.size();

  // std::cout <<" FieldPropagationHandler::DoIt(baskets) called for " << ntracks
  //          << " tracks." << std::endl;

  double *steps = td->GetDblArray(ntracks);
  for (int itr = 0; itr < ntracks; itr++) {
    // Can this loop be vectorized?
    GeantTrack &track = *tracks[itr];
    lmax = SafeLength(track, gEpsDeflection);
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
#if ! (defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  for (auto track : tracks) {
    if (track->Status() == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++;  // Find a new Counter !!! TODO
      td->fNsteps++;
      track->SetStage(stageAfterCrossing); // kPostPropagationStage);
    } else {
       // Vector treatment was not requested, so proceed with scalar
       if (!IsSameLocation(*track, td)) {
          td->fNcross++;
          td->fNsteps++;
       } else {
          track->SetStage(kGeometryStepStage);
       }
       output.AddTrack(track);
    }
  }
#else
  // If vectorized treatment was requested and the remaining population is
  // large enough, continue with vectorized treatment
  constexpr int kMinVecSize = 8; // this should be retrieved from elsewhere
  int nvect = 0;
  if (nvect < kMinVecSize) {
    for (auto track : tracks) {
      if (track->Status() == kPhysics) continue;
      if (!IsSameLocation(*track, td)) {
        td->fNcross++;
        td->fNsteps++; // Why not ?        
        track->SetStage(stageAfterCrossing); // Future: (kPostPropagationStage);
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
void FieldPropagationHandler::PropagateInVolume(GeantTrack &track, double crtstep, GeantTaskData * td)
{
// Single track propagation in a volume. The method is to be called
// only with  charged tracks in magnetic field.The method decreases the fPstepV
// fSafetyV and fSnextV with the propagated values while increasing the fStepV.
// The status and boundary flags are set according to which gets hit first:
// - physics step (bdr=0)
// - safety step (bdr=0)
// - snext step (bdr=1)

   // std::cout << "FieldPropagationHandler::PropagateInVolume called for 1 track" << std::endl;

   using ThreeVector = vecgeom::Vector3D<double>;
   bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;
   double bmag= -1.0;
   ThreeVector BfieldInitial;
   ThreeVector Position(track.X(), track.Y(), track.Z());
   FieldLookup::GetFieldValue(Position, BfieldInitial, bmag);

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
   auto fieldPropagator = GetFieldPropagator(td);
   if( !fieldPropagator || !td->fSpace4FieldProp) {
      fieldPropagator= Initialize(td);
   }
#endif

  // Reset relevant variables
  track.SetStatus(kInFlight);
  double pstep = track.GetPstep() - crtstep;
  if (pstep < 1.E-10) {
    pstep = 0;
    track.SetStatus(kPhysics);
  }
  track.SetPstep(pstep);
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

#if 0
  constexpr double inv_kilogauss = 1.0 / geant::kilogauss;
  double curvaturePlus= fabs(GeantTrack::kB2C * track.Charge() * (bmag* inv_kilogauss)) / (track.P() + 1.0e-30);  // norm for step

  const double angle= crtstep * curvaturePlus;
  constexpr double numRadiansMax= 10.0; // Too large an angle - many RK steps.  Potential change -> 2.0*PI;
  constexpr double numRadiansMin= 0.05; // Very small an angle - helix is adequate.  TBC: Use average B-field value?
      //  A track turning more than 10 radians will be treated approximately

  bool mediumAngle = ( numRadiansMin < angle ) && ( angle < numRadiansMax );
  useRungeKutta = useRungeKutta && (mediumAngle);
#endif

  // if( angle > 0.000001 ) std::cout << " ang= " << angle << std::endl;
  bool dominantBz =  std::fabs( std::fabs(BfieldInitial[2]) )
     > 1.e3 * std::max( std::fabs( BfieldInitial[0]), std::fabs(BfieldInitial[1]) );

#ifdef DEBUG_FIELD
  printf("--PropagateInVolume(Single): ");
  printf("Momentum= %9.4g (MeV) Curvature= %9.4g (1/mm)  CurvPlus= %9.4g (1/mm)  step= %f (mm)  Bmag=%8.4g KG   angle= %g\n",
         (track.P()/geant::MeV), Curvature(track)*geant::mm, curvaturePlus*geant::mm, crtstep/geant::mm,
         bmag*inv_kilogauss,  angle );
#endif

  ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
  ThreeVector PositionNew(0.,0.,0.);
  ThreeVector DirectionNew(0.,0.,0.);

// #ifndef VECCORE_CUDA
  if( useRungeKutta ) {
     fieldPropagator->DoStep(Position,    Direction,    track.Charge(), track.P(), crtstep,
                             PositionNew, DirectionNew);
#ifdef STATS_METHODS
     numRK++;
#endif
  }
  else
// #endif
  {
     constexpr double toKiloGauss= 1.0e+14; // Converts to kilogauss -- i.e. 1 / Unit::kilogauss
                                            // Must agree with values in magneticfield/inc/Units.h
     double Bz = BfieldInitial[2] * toKiloGauss;
     if ( dominantBz ) {
        // Constant field in Z-direction
        ConstBzFieldHelixStepper stepper( Bz ); //
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,    track.Charge(), track.P(), crtstep,
                                               PositionNew, DirectionNew);
#ifdef STATS_METHODS
        numHelixZ++;
#endif
     } else {
        // Geant::
        double BfieldArr[3] = { BfieldInitial.x(), BfieldInitial.y(), BfieldInitial.z() };
        ConstFieldHelixStepper stepper( BfieldArr );
        stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), crtstep,
                                               PositionNew, DirectionNew);
#ifdef STATS_METHODS
        numHelixGen++;
#endif
     }
  }
#ifdef STATS_METHODS
  unsigned long nTot;
  nTot=
     numTot++;
#ifdef PRINT_STATS
  // unsigned long nTot = numTot;
  unsigned long modbase = 10000;
  if( numTot % modbase < 1 ) {
        unsigned long rk= numRK, hZ= numHelixZ, hGen= numHelixGen;
        std::cerr << "Step statistics (field Propagation):  total= " << nTot
                  << " RK = " << rk << "  HelixGen = " << hGen << " Helix-Z = " << hZ << std::endl;
  }
  if( numTot > 10 * modbase )
     modbase = 10 * modbase;
#endif
#endif

  //  may normalize direction here  // vecCore::math::Normalize(dirnew);
  ThreeVector DirectionUnit = DirectionNew.Unit();
  double posShift = (PositionNew - Position).Mag();

  track.SetPosition(PositionNew);
  track.SetDirection(DirectionUnit);

  track.DecreaseSafety(posShift); //  Was crtstep;
  if (track.GetSafety() < 1.E-10)
    track.SetSafety(0);

#ifdef REPORT_AND_CHECK
  double origMag= Direction.Mag();
  double oldMag= DirectionNew.Mag();
  double newMag= DirectionUnit.Mag();
  Printf(" -- State after propagation in field:  Position= %f, %f, %f   Direction= %f, %f, %f  - mag original, integrated, integr-1.0, normed-1 = %10.8f %10.8f %7.2g %10.8f %7.2g",
         track.X(),  track.Y(),  track.Z(), track.Dx(),  track.Dy(),  track.Dz(),
         origMag, oldMag, oldMag-1.0, newMag, newMag-1.0 );

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
      // Geant::Print("PropagateInVolume/Single","relative difference in pos = %g", diffpos/crtstep);
      Geant::Print("PropagateInVolume/Single","difference in pos = %g (abs) %g (relative) , step= %g",
                   diffpos, diffpos/crtstep, crtstep);
  }
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PropagateInVolume(TrackVec_t &tracks,
                                                const double *stepSize,
                                                GeantTaskData *td)
{
// The Vectorized Implementation for Magnetic Field Propagation
   // std::cout << "FieldPropagationHandler::PropagateInVolume called for Many tracks" << std::endl;

  int nTracks = tracks.size();
#if 1 // VECTOR_FIELD_PROPAGATION
  using vecgeom::SOA3D;
  using vecgeom::Vector3D;
  const int Npm= 6;
  const double epsTol= 3.0e-5;

  // double yInput[8*nTracks], yOutput[8*nTracks];
  bool       succeeded[nTracks];
  int        intCharge[nTracks];
  double     fltCharge[nTracks];

  // Choice 1.  SOA3D
  PrepareBuffers(nTracks, td);

  auto  wsp = td->fSpace4FieldProp; // WorkspaceForFieldPropagation *
  SOA3D<double>& position3D  = * (wsp->fPositionInp);
  SOA3D<double>& direction3D = * (wsp->fDirectionInp);
  double        momentumMag[nTracks];

  SOA3D<double>& PositionOut  = * (wsp->fPositionOutp);
  SOA3D<double>& DirectionOut = * (wsp->fDirectionOutp);

  // Choice 2.   SOA6D
  // SOA6D<double> PositMom6D( nTracks );

  for (int itr=0; itr<nTracks; ++itr)
  {
     GeantTrack* pTrack= tracks[itr];

     intCharge[itr]= pTrack->Charge();
     fltCharge[itr]= pTrack->Charge();
     momentumMag[itr] = pTrack->P();

     // PositMom6D.push_back( pTrack->X(), pTrack->Y(), pTrack->Z(), px, py, pz );
     position3D.push_back( pTrack->X(), pTrack->Y(), pTrack->Z() );
     direction3D.push_back( pTrack->Dx(), pTrack->Dy(), pTrack->Dz() );
  }

  // Trial Helix step -- now assumes a constant uniform field
  // auto config= td->fPropagator->fConfig;
  auto fieldConfig = FieldLookup::GetFieldConfig();
  assert ( fieldConfig != nullptr);

  if( 0 ) // fieldConfig->IsFieldUniform() )
  {
     vecgeom::Vector3D<double> BfieldUniform= fieldConfig->GetUniformFieldValue();
     ConstFieldHelixStepper stepper( BfieldUniform );
     // stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), stepSize,
     //                                        PositionNew, DirectionNew);

     // std::cout << "Before Helix stepper - Position addresses: x= " << PositionOut.x() << " y= " << PositionOut.y()
     //           << " z=" << PositionOut.z() << std::endl;

     stepper.DoStepArr</*Geant::*/Double_v>( position3D.x(),  position3D.y(),  position3D.z(),
                                         direction3D.x(), direction3D.y(), direction3D.z(),
                                         intCharge,
                                         momentumMag,
                                         stepSize,
                                         PositionOut.x(),  PositionOut.y(),  PositionOut.z(),
                                         DirectionOut.x(), DirectionOut.y(), DirectionOut.z(),
                                         nTracks
        );
     // Store revised positions and location in original tracks
     for (int itr=0; itr<nTracks; ++itr)
     {
        GeantTrack& track= *tracks[itr];
        Vector3D<double> positionMove = { track.X(),   //  - PositionOut.x(itr),
                                          track.Y(),   //  - PositionOut.y(itr),
                                          track.Z() }; //  - PositionOut.z(itr) };
        positionMove -= PositionOut[itr];
        double posShift = positionMove.Mag();
        track.SetPosition(PositionOut.x(itr), PositionOut.y(itr), PositionOut.z(itr));
        track.SetDirection(DirectionOut.x(itr), DirectionOut.y(itr), DirectionOut.z(itr));

        // Check new direction
        Vector3D<double> dirOut( DirectionOut.x(itr), DirectionOut.y(itr), DirectionOut.z(itr) );
        // std::cout << "["<< itr << "] new direction = " << dirOut.x() << " " << dirOut.y() << " " << dirOut.z() << std::endl;
        assert( fabs( dirOut.Mag() - 1.0 ) < 1.0e-6 && "Out Direction is not normalised." );
        // Exact update of the safety - using true move (not distance along curve)
        track.DecreaseSafety(posShift); //  Was crtstep;
        if (track.GetSafety() < 1.E-10)
           track.SetSafety(0);
     }
  }
  else
  {
     // Prepare for Runge Kutta stepping

     bool checkPrint= false;
     if( checkPrint )
        std::cerr << " Filling fieldTrack Array (FieldPropagationHandler::PropagateInVolume (basket)" << std::endl;

     // Choice 3.   Array of FieldTrack
     FieldTrack fldTracksIn[nTracks], fldTracksOut[nTracks];
     for (int itr=0; itr<nTracks; ++itr)
     {
        GeantTrack* pTrack= tracks[itr];
        // Alternative - directly momentum vector
        double pmag= pTrack->P(), px=pTrack->Dx(), py=pTrack->Dy(), pz=pTrack->Dz();
        px *= pmag;
        py *= pmag;
        pz *= pmag;
        // Momentum3D.push( px, py, pz );

        // Choice 3. --- Load
        // yInput[itr].LoadFromTrack(*tracks[itr]);
        double trackVals[Npm] = { pTrack->X(), pTrack->Y(), pTrack->Z(), px, py, pz };
        fldTracksIn[itr].LoadFromArray( trackVals, Npm );

        if( checkPrint )
           std::cerr
              << pTrack  // tracks[itr]
              << " Pos = " << fldTracksIn[itr][0] << ", " << fldTracksIn[itr][1] << ", " << fldTracksIn[itr][2]
              << " Dir = " << fldTracksIn[itr][3] << ", " << fldTracksIn[itr][4] << ", " << fldTracksIn[itr][5]
              // EndPositionScalar
              << std::endl;
     }

     auto fieldPropagator = GetFieldPropagator(td);
     auto vectorDriver = fieldPropagator ? fieldPropagator->GetFlexibleIntegrationDriver() : nullptr;

     if( vectorDriver ) {
        // Integrate using Runge Kutta method
        vectorDriver
           ->AccurateAdvance( fldTracksIn, stepSize, fltCharge, epsTol,
                              fldTracksOut, nTracks, succeeded );

#ifdef CHECK_VS_SCALAR
        bool checkVsScalar= true;
        const char *diffBanner= "Differences (if any between vector & scalar RK method:";
        bool bannerUsed= false;
#endif
   
        // Store revised positions and location in original tracks
        for (int itr=0; itr<nTracks; ++itr)
        {
           GeantTrack& track= *tracks[itr];
           FieldTrack& fldTrackEnd= fldTracksOut[itr];
           Vector3D<double> startPosition  = { track.X(), track.Y(), track.Z() };
           Vector3D<double> startDirection = { track.Dx(), track.Dy(), track.Dz() };
           Vector3D<double> endPosition = { fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2] };
           double posShift = (startPosition-endPosition).Mag();

           // Vector3D<double> endMomentum = { fldTrackEnd[3], fldTrackEnd[4], fldTrackEnd[5] };
           const double pX= fldTrackEnd[3];
           const double pY= fldTrackEnd[4];
           const double pZ= fldTrackEnd[5];
           const double pmag_inv= 1.0 / track.P();
           // ---- Perform checks
#ifdef CHECK_VS_SCALAR           
           if( checkVsScalar )
           {
              // 1. Double check magnitude at end point
              double pMag2End = ( pX*pX + pY * pY + pZ * pZ);
              double relDiff = pMag2End * pmag_inv * pmag_inv - 1.0;
              if( std::fabs(relDiff) > geant::perMillion ) {
                 if( !bannerUsed ) { std::cerr << diffBanner << std::endl; bannerUsed= true; }                  
                 std::cerr << "Track " << itr << " has momentum magnitude difference "
                           << relDiff << "  Momentum magnitude @ end = " << std::sqrt( pMag2End )
                           << " vs. start = " << track.P() << std::endl;
              }
              assert( pMag2End > 0.0 && fabs(relDiff) < 0.01 && "ERROR in direction normal.");

              // 2. Check against 'scalar' propagation of the same track
              using ThreeVector = vecgeom::Vector3D<double>;
              ThreeVector Position(track.X(), track.Y(), track.Z());
              ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
              ThreeVector EndPositionScalar(0.,0.,0.);
              ThreeVector EndDirScalar(0.,0.,0.);
              fieldPropagator->DoStep(startPosition,     startDirection, track.Charge(), track.P(), stepSize[itr],
                                      EndPositionScalar, EndDirScalar);
              //      checking direction
              ThreeVector EndDirVector(pmag_inv * pX, pmag_inv * pY, pmag_inv * pZ);
              ThreeVector diffDir = EndDirVector - EndDirScalar;
              double diffDirMag = diffDir.Mag();
              if( diffDirMag > geant::perMillion ) {
                 if( !bannerUsed ) { std::cerr << diffBanner << std::endl; bannerUsed= true; }

                 std::cerr << "Track [" << itr << "] : direction differs " 
                           << " by " << diffDir << "  ( mag = " << diffDirMag << " ) "
                           << " Direction vector = " << EndDirVector << "  scalar = " << EndDirScalar
                           << " End position= " << endPosition                    
                           << std::endl;
              } else {
                 //   checking against magnitude of direction difference
                 ThreeVector changeDirVector = EndDirVector - startDirection;
                 ThreeVector changeDirScalar = EndDirScalar - startDirection;
                 double magChangeDirScalar = changeDirScalar.Mag();
                 const  double relativeDiffMax= 1.0e-3;
                 
                 if( diffDirMag > relativeDiffMax * magChangeDirScalar
                     && (diffDirMag         > 1.e-9)
                     && (magChangeDirScalar > 1.e-10) )
                 {
                    if( !bannerUsed ) { std::cerr << diffBanner << std::endl; bannerUsed= true; }
                    std::cerr << "Track [" << itr << "] : direction CHANGE differs"
                              << " by " << diffDir << "  ( mag = " << diffDirMag
                              << " vs |delta Dir scalar| " << magChangeDirScalar << "  ) "
                              << " DeltaDir/V = " << changeDirVector << "  scalar = " << changeDirScalar
                              << " End position= " << endPosition
                              << std::endl;
                 }
              }

              //      checking position
              ThreeVector diffPos = endPosition - EndPositionScalar;
              if( diffPos.Mag() > geant::perMillion ) {
                if( !bannerUsed ) { std::cerr << diffBanner << std::endl; bannerUsed= true; }                 
                std::cerr << "Track [" << itr << "] : position diff " << diffPos
                          << " mag= " << diffPos.Mag() << " "
                          << " Pos: vec = " << endPosition << " scalar = " << EndPositionScalar
                          << " move (vector) = " << endPosition - startPosition
                          << " its-mag= " << (endPosition - startPosition).Mag()
                          << " stepSize = " << stepSize[itr]
                          << " move (scalar) = " << EndPositionScalar - startPosition
                          << std::endl;
              }
           }
#endif /* CHECK_VS_SCALAR */
           // Update the state of this track
           track.SetPosition(fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2]);
           track.SetDirection(pmag_inv * pX, pmag_inv * pY, pmag_inv * pZ);

           // Exact update of the safety - using true move (not distance along curve)
           track.DecreaseSafety(posShift); //  Was crtstep;
           // track.DecreaseSafety(stepSize[itr]);  // Trial fix/change - in case this is sensitive
           if (track.GetSafety() < 1.E-10)
             track.SetSafety(0);
        }
     } else {
        // Geant::Error( ... );
        std::cerr << "FieldPropagationHandler: no Flexible/Vector Integration Driver found."
                  << std::endl;
     }
  }
#else
  // Placeholder - implemented just as a loop
  for (int itr=0; itr<nTracks; ++itr)
    PropagateInVolume(*tracks[itr], stepSize[itr], td);
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
