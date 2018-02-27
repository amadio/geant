#include "Geant/FieldPropagationHandler.h"

#include "Geant/FieldConfig.h"
#include "Geant/FieldLookup.h"

#include "Geant/GUFieldPropagatorPool.h"
#include "Geant/GUFieldPropagator.h"
#include "Geant/ConstBzFieldHelixStepper.h"
#include "Geant/ConstFieldHelixStepper.h"
#include "Geant/FieldTrack.h"

#include "Geant/Track.h"

#include "base/SOA3D.h"
// #include "SOA6D.h"
#include "Geant/VectorTypes.h" // Defines geant::Double_v etc
#include "Geant/SystemOfUnits.h"

#include "Geant/WorkspaceForFieldPropagation.h"
#include "Geant/FlexIntegrationDriver.h"

#include "navigation/NavigationState.h"
#include "Geant/ScalarNavInterfaceVG.h"
#include "Geant/ScalarNavInterfaceVGM.h"
#include "Geant/VectorNavInterface.h"

using Double_v = geant::Double_v;

// #define CHECK_VS_RK   1
//#define CHECK_VS_HELIX 1
#define REPORT_AND_CHECK 1

#define STATS_METHODS 1

#ifdef CHECK_VS_HELIX
#define CHECK_VS_SCALAR 1
#endif

#ifdef CHECK_VS_RK
#define CHECK_VS_SCALAR 1
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

constexpr double FieldPropagationHandler::gEpsDeflection = 1.E-2 * units::cm;

auto stageAfterCrossing = kPostPropagationStage;

#ifdef STATS_METHODS
static std::atomic<unsigned long> numRK, numHelixZ, numHelixGen, numTot;
#endif

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::FieldPropagationHandler(int threshold, Propagator *propagator, double epsTol)
    : Handler(threshold, propagator), fEpsTol(epsTol)
{
  // Default constructor
  // std::cout << " FieldPropagationHandler c-tor called:  threshold= " << threshold << std::endl;
  InitializeStats();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::InitializeStats()
{
#ifdef STATS_METHODS
  numTot      = 0;
  numRK       = 0;
  numHelixZ   = 0;
  numHelixGen = 0;
#else
  std::cout << " Field Propagation Handler: no statistics for types of steps." << std::endl;
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::~FieldPropagationHandler()
{
  // Destructor
  PrintStats();
}

//______________________________________________________________________________
GUFieldPropagator *FieldPropagationHandler::Initialize(TaskData *td)
{
  GUFieldPropagator *fieldPropagator = nullptr;
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;

  if (useRungeKutta) {
    // Initialize for the current thread -- move to Propagator::Initialize() or per thread Init method
    static GUFieldPropagatorPool *fieldPropPool = GUFieldPropagatorPool::Instance();
    assert(fieldPropPool);

    fieldPropagator = fieldPropPool->GetPropagator(td->fTid);
    assert(fieldPropagator);
    td->fFieldPropagator = fieldPropagator;
  }

  size_t basketSize    = td->fPropagator->fConfig->fNperBasket;
  td->fSpace4FieldProp = new WorkspaceForFieldPropagation(basketSize);
#endif
  // For the moment return the field Propagator.
  // Once this method is called by the framework, this can be obtained from td->fFieldPropagator
  return fieldPropagator;
}

//______________________________________________________________________________
void FieldPropagationHandler::Cleanup(TaskData *td)
{
  delete td->fSpace4FieldProp;
  td->fSpace4FieldProp = nullptr;
}

//______________________________________________________________________________
// Curvature for general field
VECCORE_ATT_HOST_DEVICE
double FieldPropagationHandler::Curvature(const Track &track) const
{
  using ThreeVector_d            = vecgeom::Vector3D<double>;
  constexpr double tiny          = 1.E-30;
  constexpr double inv_kilogauss = 1.0 / units::kilogauss;
  ThreeVector_d MagFld;
  double bmag = 0.0;

  ThreeVector_d Position(track.X(), track.Y(), track.Z());
  FieldLookup::GetFieldValue(Position, MagFld, bmag); // , td);
  bmag *= inv_kilogauss;

  //  Calculate transverse momentum 'Pt' for field 'B'
  //
  ThreeVector_d Momentum(track.Px(), track.Py(), track.Pz());
  ThreeVector_d PtransB; //  Transverse wrt direction of B
  double ratioOverFld        = 0.0;
  if (bmag > 0) ratioOverFld = Momentum.Dot(MagFld) / (bmag * bmag);
  PtransB                    = Momentum - ratioOverFld * MagFld;
  double Pt_mag              = PtransB.Mag();

  return fabs(Track::kB2C * track.Charge() * bmag / (Pt_mag + tiny));
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::DoIt(Track *track, Basket &output, TaskData *td)
{
  // Scalar geometry length computation. The track is moved into the output basket.
  // Step selection
  double step, lmax;

  // std::cout <<" FieldPropagationHandler::DoIt(*track) called for 1 ptrTrack." << std::endl;

  // We use the track sagitta to estimate the "bending" error,
  // i.e. what is the propagated length for which the track deviation in
  // magnetic field with respect to straight propagation is less than epsilon.
  // Take the maximum between the safety and the "bending" safety
  lmax = SafeLength(*track, gEpsDeflection);
  lmax = vecCore::math::Max<double>(lmax, track->GetSafety());
  // Select step to propagate as the minimum among the "safe" step and:
  // the straight distance to boundary (if frombdr=1) or the proposed  physics
  // step (frombdr=0)
  step = (track->Boundary()) ? vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track->GetSnext(), 1.E-4))
                             : vecCore::math::Min<double>(lmax, track->GetPstep());
  // Propagate in magnetic field
  PropagateInVolume(*track, step, td);
  // Update number of partial steps propagated in field
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
void FieldPropagationHandler::DoIt(Basket &input, Basket &output, TaskData *td)
{
  // Vector geometry length computation. The tracks are moved into the output basket.
  TrackVec_t &tracks = input.Tracks();
  double lmax;
  // const double gEpsDeflection = 1.E-2 * geant::units::cm;  // Units!

  // using minD= vecCore::math::Min<double>;
  // using maxD= vecCore::math::Max<double>;

  int ntracks = tracks.size();

  // std::cout <<" FieldPropagationHandler::DoIt(baskets) called for " << ntracks
  //          << " tracks." << std::endl;

  double *steps = td->GetDblArray(ntracks);
  for (int itr = 0; itr < ntracks; itr++) {
    // Can this loop be vectorized?
    Track &track = *tracks[itr];
    lmax         = SafeLength(track, gEpsDeflection);
    lmax         = vecCore::math::Max<double>(lmax, track.GetSafety());
    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = (track.Boundary())
                     ? vecCore::math::Min<double>(lmax, vecCore::math::Max<double>(track.GetSnext(), 1.E-4))
                     : vecCore::math::Min<double>(lmax, track.GetPstep());
  }
  // Propagate the vector of tracks
  PropagateInVolume(input.Tracks(), steps, td);

  // Update number of partial steps propagated in field
  td->fNmag += ntracks;

  // Update time of flight and number of interaction lengths.
  // Check also if it makes sense to call the vector interfaces
  
#if !(defined(VECTORIZED_GEOMERY) && defined(VECTORIZED_SAMELOC))
  for (auto track : tracks) {
    if (track->Status() == kPhysics) {
      // Update number of steps to physics and total number of steps
      td->fNphys++; // Find a new Counter !!! TODO
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
    }
    output.AddTrack(track);
  }
#else
  // If vectorized treatment was requested and the remaining population is
  // large enough, continue with vectorized treatment
  constexpr int kMinVecSize = 8; // this should be retrieved from elsewhere
  int nvect                 = 0;
  if (nvect < kMinVecSize) {
    for (auto track : tracks) {
      if (track->Status() == kPhysics) {
        output.AddTrack(track);
        continue;
      }
      if (!IsSameLocation(*track, td)) {
        td->fNcross++;
        td->fNsteps++;                       // Why not ?
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
  TrackGeo_v &track_geo = *td.fGeoTrack;
  for (auto track : tracks) {
    if (track.Status() != kPhysics && (track.GetSafety() < 1.E-10 || track.GetSnext() < 1.E-10))
      track_geo.AddTrack(*track);
  }
  bool *same                = td->GetBoolArray(nvect);
  NavigationState *tmpstate = td->GetPath();
  VectorNavInterface::NavIsSameLocation(nvect, track_geo.fXposV, track_geo.fYposV, track_geo.fZposV, track_geo.fXdirV,
                                        track_geo.fYdirV, track_geo.fZdirV, (const VolumePath_t **)fPathV, fNextpathV,
                                        same, tmpstate);
  track_geo.UpdateOriginalTracks();
  for (itr = 0; itr < nsel; itr++) {
    Track *track = track_geo.fOriginalV[itr];
    if (!same[itr]) {
      td->fNcross++;
      td->fNsteps++;
      track->SetBoundary(true);
      track->SetStatus(kBoundary);
      if (track->NextPath()->IsOutside()) track->SetStatus(kExitingSetup);
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
void FieldPropagationHandler::PropagateInVolume(Track &track, double crtstep, TaskData *td)
{
  // Single track propagation in a volume. The method is to be called
  // only with  charged tracks in magnetic field.The method decreases the fPstepV
  // fSafetyV and fSnextV with the propagated values while increasing the fStepV.
  // The status and boundary flags are set according to which gets hit first:
  // - physics step (bdr=0)
  // - safety step (bdr=0)
  // - snext step (bdr=1)

  // std::cout << "FieldPropagationHandler::PropagateInVolume called for 1 track" << std::endl;

  using ThreeVector            = vecgeom::Vector3D<double>;
  constexpr double toKiloGauss = 1.0 / units::kilogauss; // Converts to kilogauss
  bool useRungeKutta           = td->fPropagator->fConfig->fUseRungeKutta;
  double bmag                  = -1.0;
  ThreeVector BfieldInitial;
  ThreeVector Position(track.X(), track.Y(), track.Z());
  FieldLookup::GetFieldValue(Position, BfieldInitial, bmag);

#ifndef VECCORE_CUDA_DEVICE_COMPILATION
  auto fieldPropagator = GetFieldPropagator(td);
  if (!fieldPropagator || !td->fSpace4FieldProp) {
    fieldPropagator = Initialize(td);
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
    if (track.Boundary()) track.SetStatus(kBoundary);
  }
  track.SetSnext(snext);
  track.IncreaseStep(crtstep);

#if 0
  constexpr double numRadiansMax= 10.0; // Too large an angle - many RK steps.  Potential change -> 2.0*PI;
  constexpr double numRadiansMin= 0.05; // Very small an angle - helix is adequate.  TBC: Use average B-field value?
      //  A track turning more than 10 radians will be treated approximately

  bool mediumAngle = ( numRadiansMin < angle ) && ( angle < numRadiansMax );
  useRungeKutta = useRungeKutta && (mediumAngle);
#endif

// if( angle > 0.000001 ) std::cout << " ang= " << angle << std::endl;
#ifdef DEBUG_FIELD
  double curvaturePlus= fabs(Track::kB2C * track.Charge() * (bmag * toKiloGauss)) / (track.P() + 1.0e-30);  // norm for step

  const double angle = crtstep * curvaturePlus;
  printf("--PropagateInVolume(Single): ");
  printf("Momentum= %9.4g (MeV) Curvature= %9.4g (1/mm)  CurvaturePlus = %9.4g (1/mm)  step= %f (mm)  Bmag=%8.4g KG   angle= %g\n",
         (track.P() / units::MeV), Curvature(track) * units::mm, curvaturePlus * units::mm,
         crtstep / units::mm, bmag * toKiloGauss, angle);
#endif

  ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
  ThreeVector PositionNew(0., 0., 0.);
  ThreeVector DirectionNew(0., 0., 0.);

  if (useRungeKutta) {
    fieldPropagator->DoStep(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew, DirectionNew);
// cross check
#ifdef DEBUG_FIELD
    ThreeVector PositionNewCheck(0., 0., 0.);
    ThreeVector DirectionNewCheck(0., 0., 0.);
    ConstFieldHelixStepper stepper(BfieldInitial * toKiloGauss);
    stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNewCheck, DirectionNewCheck);

    double posShift = (PositionNew - PositionNewCheck).Mag();
    double dirShift = (DirectionNew - DirectionNewCheck).Mag();
    if (posShift > 1.e-6 || dirShift > 1.e-6) {
      std::cout << "*** position/direction shift RK vs. HelixStepper :" << posShift << " / " << dirShift << "\n";
    }
#endif

#ifdef STATS_METHODS
    numRK++;
#endif
  } else {
    // Must agree with values in magneticfield/inc/Units.h
    double Bz             = BfieldInitial[2] * toKiloGauss;
    const bool dominantBz = false; // std::fabs(std::fabs(BfieldInitial[2])) >
    // 1.e3 * std::max(std::fabs(BfieldInitial[0]), std::fabs(BfieldInitial[1]));
    if (dominantBz) {
      // Constant field in Z-direction
      ConstBzFieldHelixStepper stepper(Bz); //
      stepper.DoStep<ThreeVector, double, int>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew,
                                               DirectionNew);
#ifdef STATS_METHODS
      numHelixZ++;
#endif
    } else {
      // geant::
      double BfieldArr[3] = {BfieldInitial.x(), BfieldInitial.y(), BfieldInitial.z()};
      ConstFieldHelixStepper stepper(BfieldArr);
      stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew, DirectionNew);
#ifdef STATS_METHODS
      numHelixGen++;
#endif
    }
  }

#ifdef STATS_METHODS
  unsigned long modbase = 10000;
  if (numTot % modbase < 1) {
    PrintStats();
    if (numTot > 10 * modbase) modbase = 10 * modbase;
  }
#endif

  //  may normalize direction here  // vecCore::math::Normalize(dirnew);
  ThreeVector DirectionUnit = DirectionNew.Unit();
  double posShift           = (PositionNew - Position).Mag();

  track.SetPosition(PositionNew);
  track.SetDirection(DirectionUnit);

  track.DecreaseSafety(posShift); //  Was crtstep;
  if (track.GetSafety() < 1.E-10) track.SetSafety(0);

#ifdef REPORT_AND_CHECK
  /*****
  double origMag= Direction.Mag();
  double oldMag= DirectionNew.Mag();
  double newMag= DirectionUnit.Mag();
  Printf(" -- State after propagation in field:  event %d track %p Position= %f, %f, %f   Direction= %f, %f, %f  - mag
  original, integrated, integr-1.0, normed-1 = %7.5f %7.5f %4.2g %7.5f %4.2g", track.Event(), (void*) &track, track.X(),
  track.Y(),  track.Z(), track.Dx(),  track.Dy(),  track.Dz(), origMag, oldMag, oldMag-1.0, newMag, newMag-1.0 );
  ****/

  // int propagationType = 5;
  /*  const char* Msg[4]= { "After propagation in field - type Unknown(ERROR) ",
                        "After propagation in field - with RK           ",
                        "After propagation in field - with Helix-Bz     ",
                        "After propagation in field - with Helix-General" };  */

  CheckTrack(track, "End of Propagate-In-Volume", 1.0e-5); // Msg[propagationType] );
#endif

#if 0
  ThreeVector SimplePosition = Position + crtstep * Direction;
  // double diffpos2 = (PositionNew - Position).Mag2();
  double diffpos2 = (PositionNew - SimplePosition).Mag2();
  //   -- if (vecCore::math::Sqrt(diffpos)>0.01*crtstep) {
  const double drift= 0.01*crtstep;
  if ( diffpos2>drift*drift ){
      double diffpos= vecCore::math::Sqrt(diffpos2);
      // geant::Print("PropagateInVolume/Single","relative difference in pos = %g", diffpos/crtstep);
      geant::Print("PropagateInVolume/Single","difference in pos = %g (abs) %g (relative) , step= %g",
                   diffpos, diffpos/crtstep, crtstep);
  }
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PropagateInVolume(TrackVec_t &tracks, const double *stepSize, TaskData *td)
{
  // The Vectorized Implementation for Magnetic Field Propagation
  // std::cout << "FieldPropagationHandler::PropagateInVolume called for Many tracks" << std::endl;

  using ThreeVector            = vecgeom::Vector3D<double>;
  constexpr double toKiloGauss = 1.0 / units::kilogauss; // Converts to kilogauss
  int nTracks                  = tracks.size();
#if 1 // VECTOR_FIELD_PROPAGATION
  using vecgeom::SOA3D;
  using vecgeom::Vector3D;
  const int Npm = 6;

  // double yInput[8*nTracks], yOutput[8*nTracks];
  bool succeeded[nTracks];
  // int        intCharge[nTracks];

  // Choice 1.  SOA3D
  PrepareBuffers(nTracks, td);

  auto wsp                   = td->fSpace4FieldProp; // WorkspaceForFieldPropagation *
  double *fltCharge          = wsp->fChargeInp;
  double *momentumMag        = wsp->fMomentumInp;
  double *steps              = wsp->fStepsInp;
  SOA3D<double> &position3D  = *(wsp->fPositionInp);
  SOA3D<double> &direction3D = *(wsp->fDirectionInp);

  SOA3D<double> &PositionOut  = *(wsp->fPositionOutp);
  SOA3D<double> &DirectionOut = *(wsp->fDirectionOutp);

  // Choice 2.   SOA6D
  // SOA6D<double> PositMom6D( nTracks );

  for (int itr = 0; itr < nTracks; ++itr) {
    Track *pTrack = tracks[itr];

    // intCharge[itr]= pTrack->Charge();
    fltCharge[itr]   = pTrack->Charge();
    momentumMag[itr] = pTrack->P();
    steps[itr]       = stepSize[itr];

    // PositMom6D.push_back( pTrack->X(), pTrack->Y(), pTrack->Z(), px, py, pz );
    position3D.push_back(pTrack->X(), pTrack->Y(), pTrack->Z());
    direction3D.push_back(pTrack->Dx(), pTrack->Dy(), pTrack->Dz());
  }

  // Trial Helix step -- now assumes a constant uniform field
  // auto config= td->fPropagator->fConfig;
  auto fieldConfig = FieldLookup::GetFieldConfig();
  assert(fieldConfig != nullptr);

  if (0 /*fieldConfig->IsFieldUniform()*/) {
    vecgeom::Vector3D<double> BfieldUniform = fieldConfig->GetUniformFieldValue();
    ConstFieldHelixStepper stepper(BfieldUniform * toKiloGauss);
    // stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), stepSize,
    //                                        PositionNew, DirectionNew);

    // std::cout << "Before Helix stepper - Position addresses: x= " << PositionOut.x() << " y= " << PositionOut.y()
    //           << " z=" << PositionOut.z() << std::endl;

    stepper.DoStepArr<Double_v>(position3D.x(), position3D.y(), position3D.z(), direction3D.x(), direction3D.y(),
                                direction3D.z(), fltCharge, momentumMag, steps, PositionOut.x(), PositionOut.y(),
                                PositionOut.z(), DirectionOut.x(), DirectionOut.y(), DirectionOut.z(), nTracks);
    // Store revised positions and location in original tracks
    for (int itr = 0; itr < nTracks; ++itr) {
      Track &track = *tracks[itr];

      // Crosscheck against helix stepper
      ThreeVector Position(track.X(), track.Y(), track.Z());
      ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
#ifdef DEBUG_FIELD
      ThreeVector PositionNew(0., 0., 0.);
      ThreeVector DirectionNew(0., 0., 0.);

      stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), stepSize[itr], PositionNew, DirectionNew);
      double posDiff = (PositionNew - PositionOut[itr]).Mag();
      double dirDiff = (DirectionNew - DirectionOut[itr]).Mag();
      if (posDiff > 1.e-6 || dirDiff > 1.e-6) {
        std::cout << "*** position/direction shift HelixStepper scalar vs. vector :" << posDiff << " / " << dirDiff << "\n";
        stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), steps[itr], PositionNew, DirectionNew);
      }
#endif
      Vector3D<double> positionMove = {track.X(),  //  - PositionOut.x(itr),
                                       track.Y(),  //  - PositionOut.y(itr),
                                       track.Z()}; //  - PositionOut.z(itr) };
      positionMove -= PositionOut[itr];
      double posShift = positionMove.Mag();
      track.SetPosition(PositionOut.x(itr), PositionOut.y(itr), PositionOut.z(itr));
      track.SetDirection(DirectionOut.x(itr), DirectionOut.y(itr), DirectionOut.z(itr));
      // Update status, step and safety
      track.SetStatus(kInFlight);
      double pstep = track.GetPstep() - stepSize[itr];
      if (pstep < 1.E-10) {
        pstep = 0;
        track.SetStatus(kPhysics);
      }
      double snext = track.GetSnext() - stepSize[itr];
      if (snext < 1.E-10) {
        snext = 0;
        if (track.Boundary()) track.SetStatus(kBoundary);
      }
      track.SetSnext(snext);
      track.IncreaseStep(stepSize[itr]);

      // Check new direction
      Vector3D<double> dirOut(DirectionOut.x(itr), DirectionOut.y(itr), DirectionOut.z(itr));
      // std::cout << "["<< itr << "] new direction = " << dirOut.x() << " " << dirOut.y() << " " << dirOut.z() <<
      // std::endl;
      assert(fabs(dirOut.Mag() - 1.0) < 1.0e-6 && "Out Direction is not normalised.");
      // Exact update of the safety - using true move (not distance along curve)
      track.DecreaseSafety(posShift); //  Was crtstep;
      if (track.GetSafety() < 1.E-10) track.SetSafety(0);
    }
  } else {
    // Prepare for Runge Kutta stepping

    bool checkPrint = false;
    if (checkPrint)
      std::cerr << " Filling fieldTrack Array (FieldPropagationHandler::PropagateInVolume (basket)" << std::endl;

    // Choice 3.   Array of FieldTrack
    FieldTrack fldTracksIn[nTracks], fldTracksOut[nTracks];
    for (int itr = 0; itr < nTracks; ++itr) {
      Track *pTrack = tracks[itr];
      // Alternative - directly momentum vector
      double pmag = pTrack->P(), px = pTrack->Dx(), py = pTrack->Dy(), pz = pTrack->Dz();
      px *= pmag;
      py *= pmag;
      pz *= pmag;
      // Momentum3D.push( px, py, pz );

      // Choice 3. --- Load
      // yInput[itr].LoadFromTrack(*tracks[itr]);
      double trackVals[Npm] = {pTrack->X(), pTrack->Y(), pTrack->Z(), px, py, pz};
      fldTracksIn[itr].LoadFromArray(trackVals, Npm);

      if (checkPrint)
        std::cerr << pTrack // tracks[itr]
                  << " Pos = " << fldTracksIn[itr][0] << ", " << fldTracksIn[itr][1] << ", " << fldTracksIn[itr][2]
                  << " Dir = " << fldTracksIn[itr][3] << ", " << fldTracksIn[itr][4] << ", " << fldTracksIn[itr][5]
                  // EndPositionScalar
                  << std::endl;
    }

    auto fieldPropagator = GetFieldPropagator(td);
    auto vectorDriver    = fieldPropagator ? fieldPropagator->GetFlexibleIntegrationDriver() : nullptr;

    if (vectorDriver) {
      // Integrate using Runge Kutta method
      vectorDriver->AccurateAdvance(fldTracksIn, steps, fltCharge, fEpsTol, fldTracksOut, nTracks, succeeded);

#ifdef CHECK_VS_SCALAR
      bool checkVsScalar     = true;
      const char *diffBanner = "Differences (if any between vector & scalar RK method:";
      bool bannerUsed        = false;
#endif

      // Store revised positions and location in original tracks
      for (int itr = 0; itr < nTracks; ++itr) {
        Track &track                    = *tracks[itr];
        FieldTrack &fldTrackEnd         = fldTracksOut[itr];
        Vector3D<double> startPosition  = {track.X(), track.Y(), track.Z()};
        Vector3D<double> startDirection = {track.Dx(), track.Dy(), track.Dz()};
        Vector3D<double> endPosition    = {fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2]};
        double posShift                 = (startPosition - endPosition).Mag();

        // Vector3D<double> endMomentum = { fldTrackEnd[3], fldTrackEnd[4], fldTrackEnd[5] };
        const double pX       = fldTrackEnd[3];
        const double pY       = fldTrackEnd[4];
        const double pZ       = fldTrackEnd[5];
        const double pmag_inv = 1.0 / track.P();

        // ---- Perform checks
        using ThreeVector = vecgeom::Vector3D<double>;
        ThreeVector endDirVector(pmag_inv * pX, pmag_inv * pY, pmag_inv * pZ);
        endDirVector.Normalize();

        ThreeVector endPositionScalar(0., 0., 0.), endDirScalar(0., 0., 0.);
        fieldPropagator->DoStep(startPosition, startDirection, track.Charge(), track.P(), stepSize[itr], endPositionScalar,
                                  endDirScalar);
        double posErr = (endPositionScalar - endPosition).Mag();
        double dirErr = (endDirScalar - endDirVector).Mag();
        if (posErr > 1.e-3 * posShift || dirErr > 1.e-6) {
          std::cout << "*** position/direction shift scalar RK vs. vector RK :" << posErr << " / " << dirErr << "\n";
        }

#ifdef CHECK_VS_SCALAR
        if (checkVsScalar) {
          // 1. Double check magnitude at end point
          double pMag2End = (pX * pX + pY * pY + pZ * pZ);
          double relDiff  = pMag2End * pmag_inv * pmag_inv - 1.0;
          if (std::fabs(relDiff) > units::perMillion) {
            if (!bannerUsed) {
              std::cerr << diffBanner << std::endl;
              bannerUsed = true;
            }
            std::cerr << "Track " << itr << " has momentum magnitude difference " << relDiff
                      << "  Momentum magnitude @ end = " << std::sqrt(pMag2End) << " vs. start = " << track.P()
                      << std::endl;
          }
          assert(pMag2End > 0.0 && fabs(relDiff) < 0.01 && "ERROR in direction normal.");

          // 2. Check against 'scalar' propagation of the same track
          using ThreeVector = vecgeom::Vector3D<double>;

          ThreeVector EndPositionScalar(0., 0., 0.), EndDirScalar(0., 0., 0.);

#ifdef CHECK_VS_HELIX
          vecgeom::Vector3D<double> BfieldInitial;
          double bmag, charge = track.Charge();
          FieldLookup::GetFieldValue(startPosition, BfieldInitial, bmag);
          ConstFieldHelixStepper stepper(BfieldInitial);
          stepper.DoStep<double>(startPosition, startDirection, charge, track.P(), stepSize[itr], EndPositionScalar,
                                 EndDirScalar);
#else
          fieldPropagator->DoStep(startPosition, startDirection, charge, track.P(), stepSize[itr], EndPositionScalar,
                                  EndDirScalar);
#endif
          //      checking direction
          ThreeVector EndDirVector(pmag_inv * pX, pmag_inv * pY, pmag_inv * pZ);
          ThreeVector diffDir     = EndDirVector - EndDirScalar;
          double diffDirMag       = diffDir.Mag();
          const double maxDiffMom = 1.5 * fEpsTol; // 10.0 * units::perMillion;
          if (diffDirMag > maxDiffMom) {
            if (!bannerUsed) {
              std::cerr << diffBanner << std::endl;
              bannerUsed = true;
            }
            // const Track* pTrack= tracks[itr];
            double curv = Curvature(track); // (*tracks[itr] ) ;
            std::cerr << "Track [" << itr << "] : direction differs "
                      << " by " << diffDir << "  ( mag = " << diffDirMag << " ) "
                      << " Direction vector = " << EndDirVector << "  scalar = " << EndDirScalar
                      << " stepSize = " << stepSize[itr] << " curv = " << curv << " End position= " << endPosition
                      << std::endl;
          } else {
            //   checking against magnitude of direction difference
            ThreeVector changeDirVector  = EndDirVector - startDirection;
            ThreeVector changeDirScalar  = EndDirScalar - startDirection;
            double magChangeDirScalar    = changeDirScalar.Mag();
            const double relativeDiffMax = 1.0e-3;

            if (diffDirMag > relativeDiffMax * magChangeDirScalar && (diffDirMag > 1.e-9) &&
                (magChangeDirScalar > 1.e-10)) {
              if (!bannerUsed) {
                std::cerr << diffBanner << std::endl;
                bannerUsed = true;
              }
              std::cerr << "Track [" << itr << "] : direction CHANGE differs"
                        << " by " << diffDir << "  ( mag = " << diffDirMag << " vs |delta Dir scalar| "
                        << magChangeDirScalar << "  ) "
                        << " DeltaDir/V = " << changeDirVector << "  scalar = " << changeDirScalar
                        << " End position= " << endPosition << std::endl;
            }
          }

          //      checking position
          ThreeVector diffPos     = endPosition - EndPositionScalar;
          const double maxDiffPos = 1.5 * fEpsTol; // * distanceAlongPath
          if (diffPos.Mag() > maxDiffPos) {
            if (!bannerUsed) {
              std::cerr << diffBanner << std::endl;
              bannerUsed = true;
            }
            std::cerr << "Track [" << itr << "] : position diff " << diffPos << " mag= " << diffPos.Mag() << " "
                      << " Pos: vec = " << endPosition << " scalar = " << EndPositionScalar
                      << " move (vector) = " << endPosition - startPosition
                      << " its-mag= " << (endPosition - startPosition).Mag() << " stepSize = " << stepSize[itr]
                      << " move (scalar) = " << EndPositionScalar - startPosition << std::endl;
          }
        }
#endif /* CHECK_VS_SCALAR */
        // Update the state of this track
        track.SetPosition(fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2]);
        track.SetDirection(pmag_inv * pX, pmag_inv * pY, pmag_inv * pZ);
        track.SetStatus(kInFlight);
        double pstep = track.GetPstep() - stepSize[itr];
        if (pstep < 1.E-10) {
          pstep = 0;
          track.SetStatus(kPhysics);
        }
        double snext = track.GetSnext() - stepSize[itr];
        if (snext < 1.E-10) {
          snext = 0;
          if (track.Boundary()) track.SetStatus(kBoundary);
        }
        track.SetSnext(snext);
        track.IncreaseStep(stepSize[itr]);

        // Exact update of the safety - using true move (not distance along curve)
        track.DecreaseSafety(posShift); //  Was crtstep;
        // track.DecreaseSafety(stepSize[itr]);  // Trial fix/change - in case this is sensitive
        if (track.GetSafety() < 1.E-10) track.SetSafety(0);
      }
    } else {
      // geant::Error( ... );
      std::cerr << "FieldPropagationHandler: no Flexible/Vector Integration Driver found." << std::endl;
    }
  }
#else
  // Placeholder - implemented just as a loop
  for (int itr = 0; itr < nTracks; ++itr)
    PropagateInVolume(*tracks[itr], stepSize[itr], td);
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
bool FieldPropagationHandler::IsSameLocation(Track &track, TaskData *td)
{
  // Query geometry if the location has changed for a track
  // Returns number of tracks crossing the boundary (0 or 1)

  if (track.GetSafety() > 1.E-10 && track.GetSnext() > 1.E-10) {
    // Track stays in the same volume
    track.SetBoundary(false);
    return true;
  }

  // Track may have crossed, check it
  bool same;
  vecgeom::NavigationState *tmpstate = td->GetPath();
  ScalarNavInterfaceVGM::NavIsSameLocation(track, same, tmpstate);
  if (same) {
    track.SetBoundary(false);
    return true;
  }

  track.SetBoundary(true);
  track.SetStatus(kBoundary);
  if (track.NextPath()->IsOutside()) track.SetStatus(kExitingSetup);
  if (track.GetStep() < 1.E-8) td->fNsmall++;
  return false;
}

#define IsNan(x) (!(x > 0 || x <= 0.0))
//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::CheckTrack(Track &track, const char *msg, double epsilon) const
{
  // Ensure that values are 'sensible' - else print msg and track
  if (epsilon <= 0.0 || epsilon > 0.01) {
    epsilon = 1.e-6;
  }

  double x = track.X(), y = track.Y(), z = track.Z();
  bool badPosition       = IsNan(x) || IsNan(y) || IsNan(z);
  const double maxRadius = 10000.0; // Should be a property of the geometry
  const double maxRadXY  = 5000.0;  // Should be a property of the geometry

  // const double maxUnitDev =  1.0e-4;  // Deviation from unit of the norm of the direction
  double radiusXy2 = x * x + y * y;
  double radius2   = radiusXy2 + z * z;
  badPosition      = badPosition || (radiusXy2 > maxRadXY * maxRadXY) || (radius2 > maxRadius * maxRadius);

  const double maxUnitDev = epsilon; // Use epsilon for max deviation of direction norm from 1.0

  double dx = track.Dx(), dy = track.Dy(), dz = track.Dz();
  double dirNorm2   = dx * dx + dy * dy + dz * dz;
  bool badDirection = std::fabs(dirNorm2 - 1.0) > maxUnitDev;
  if (badPosition || badDirection) {
    static const char *errMsg[4] = {" All ok - No error. ",
                                    " Bad position.",                 // [1]
                                    " Bad direction.",                // [2]
                                    " Bad direction and position. "}; // [3]
    int iM                       = 0;
    if (badPosition) {
      iM++;
    }
    if (badDirection) {
      iM += 2;
    }
    // if( badDirection ) {
    //   Printf( " Norm^2 direction= %f ,  Norm -1 = %g", dirNorm2, sqrt(dirNorm2)-1.0 );
    // }
    Printf("ERROR> Problem with track %p . Issue: %s. Info message: %s -- Mag^2(dir)= %9.6f Norm-1= %g", (void *)&track,
           errMsg[iM], msg, dirNorm2, sqrt(dirNorm2) - 1.0);
    track.Print(msg);
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PrintStats()
{
#ifdef STATS_METHODS
  unsigned long nTot = numTot++;
  unsigned long rk = numRK, hZ = numHelixZ, hGen = numHelixGen;
  std::cerr << "Step statistics (field Propagation):  total= " << nTot << " RK = " << rk << "  HelixGen = " << hGen
            << " Helix-Z = " << hZ << std::endl;
#endif
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
