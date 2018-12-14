#include "Geant/FieldPropagationHandler.h"

#include "Geant/FieldConfig.h"
#include "Geant/FieldLookup.h"

#include "Geant/GUFieldPropagatorPool.h"
#include "Geant/GUFieldPropagator.h"
#include "Geant/ConstBzFieldHelixStepper.h"
#include "Geant/ConstFieldHelixStepper.h"
#include "Geant/FieldTrack.h"

#include "Geant/Track.h"

#include <sstream>
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

// #define REPORT_AND_CHECK 1

// #define STATS_METHODS 1
// #define DEBUG_FIELD   1

#ifdef CHECK_VS_HELIX
#define CHECK_VS_SCALAR 1
#endif

#ifdef CHECK_VS_RK
#define CHECK_VS_SCALAR 1
#endif

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

const double FieldPropagationHandler::gEpsDeflection = 1.E-2 * units::cm;

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
  SetName("FieldPropagation");
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
// std::cout << " Field Propagation Handler: no statistics for types of steps." << std::endl;
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
FieldPropagationHandler::~FieldPropagationHandler()
{
  // Destructor
  std::cerr << "Statistics from FieldPropagation destructor." << std::endl;
  PrintStats();
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
  MagFld *= inv_kilogauss;

  //  Calculate transverse momentum 'Pt' for field 'B'
  //
  ThreeVector_d Momentum(track.Px(), track.Py(), track.Pz());
  ThreeVector_d PtransB; //  Transverse wrt direction of B
  double ratioOverFld = 0.0;
  if (bmag > 0) ratioOverFld = Momentum.Dot(MagFld) / (bmag * bmag);
  PtransB       = Momentum - ratioOverFld * MagFld;
  double Pt_mag = PtransB.Mag();

  return fabs(Track::kB2C * track.Charge() * bmag / (Pt_mag + tiny));
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::DoIt(Track *track, Basket &output, TaskData *td)
{
  // Scalar geometry length computation. The track is moved into the output basket.
  using vecCore::math::Max;
  using vecCore::math::Min;
  constexpr double step_push = 1.e-4;
  // The minimum step is step_push in case the physics step limit is not smaller
  double step_min = Min(track->GetPstep(), step_push);
  // The track snext value is already the minimum between geometry and physics
  double step_geom_phys = Max(step_min, track->GetSnext());
  // Field step limit. We use the track sagitta to estimate the "bending" error,
  // i.e. what is the propagated length for which the track deviation in
  // magnetic field with respect to straight propagation is less than epsilon.
  double step_field = Max(SafeLength(*track, gEpsDeflection), track->GetSafety());

  double step = Min(step_geom_phys, step_field);
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
  using vecCore::math::Max;
  using vecCore::math::Min;
  constexpr double step_push = 1.e-4;
  TrackVec_t &tracks         = input.Tracks();
  int ntracks                = tracks.size();

  double *steps = td->GetDblArray(ntracks);
  for (int itr = 0; itr < ntracks; itr++) {
    // Can this loop be vectorized?
    Track &track = *tracks[itr];
    // The minimum step is step_push in case the physics step limit is not smaller
    double step_min = Min(track.GetPstep(), step_push);
    // The track snext value is already the minimum between geometry and physics
    double step_geom_phys = Max(step_min, track.GetSnext());
    // Field step limit. We use the track sagitta to estimate the "bending" error,
    // i.e. what is the propagated length for which the track deviation in
    // magnetic field with respect to straight propagation is less than epsilon.
    double step_field = Max(SafeLength(track, gEpsDeflection), track.GetSafety());

    // Select step to propagate as the minimum among the "safe" step and:
    // the straight distance to boundary (if fboundary=1) or the proposed  physics
    // step (fboundary=0)
    steps[itr] = Min(step_geom_phys, step_field);
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
      // if (track->GetStep() < 1.E-8) td->fNsmall++;
    } else {
      track->SetBoundary(false);
      track->SetStage(kGeometryStepStage);
    }
    output.AddTrack(track);
  }

#endif
}

void FieldPropagationHandler::VectorDispatchOverhead(TrackVec_t &tracks, TaskData *td)
{
  // This function just reproduces the overhead for gather/scatter in vector mode
  using vecgeom::SOA3D;
  using vecgeom::Vector3D;
  const int nTracks = tracks.size();

  PrepareBuffers(nTracks, td);

  auto wsp                   = td->fSpace4FieldProp; // WorkspaceForFieldPropagation *
  double *fltCharge          = wsp->fChargeInp;
  double *momentumMag        = wsp->fMomentumInp;
  double *steps              = wsp->fStepsInp;
  SOA3D<double> &position3D  = *(wsp->fPositionInp);
  SOA3D<double> &direction3D = *(wsp->fDirectionInp);

  SOA3D<double> &PositionOut  = *(wsp->fPositionOutp);
  SOA3D<double> &DirectionOut = *(wsp->fDirectionOutp);

  for (int itr = 0; itr < nTracks; ++itr) {
    Track *pTrack = tracks[itr];
    // gather
    fltCharge[itr]   = pTrack->Charge();
    momentumMag[itr] = pTrack->P();
    steps[itr]       = 0;

    position3D.push_back(pTrack->X(), pTrack->Y(), pTrack->Z());
    direction3D.push_back(pTrack->Dx(), pTrack->Dy(), pTrack->Dz());
    // scatter
    PositionOut.push_back(pTrack->X(), pTrack->Y(), pTrack->Z());
    DirectionOut.push_back(pTrack->Dx(), pTrack->Dy(), pTrack->Dz());
  }
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

  bool useRungeKutta = td->fPropagator->fConfig->fUseRungeKutta;
  auto fieldConfig   = FieldLookup::GetFieldConfig();
  double bmag        = -1.0;

  ThreeVector BfieldInitial;
  ThreeVector Position(track.X(), track.Y(), track.Z());
  FieldLookup::GetFieldValue(Position, BfieldInitial, bmag);

  auto fieldPropagator = td->fFieldPropagator;

#if DEBUG_FIELD
  bool verboseDiff = true; // If false, print just one line.  Else more details.
  bool epsilonRK   = td->fPropagator->fConfig->fEpsilonRK;
  double curvaturePlus =
      fabs(Track::kB2C * track.Charge() * (bmag * toKiloGauss)) / (track.P() + 1.0e-30); // norm for step
  const double angle = crtstep * curvaturePlus;
#endif

#ifdef PRINT_STEP_SINGLE
  Print("--PropagateInVolume(Single): ",
        "Momentum= %9.4g (MeV) Curvature= %9.4g (1/mm)  CurvPlus= %9.4g (1/mm)  step= "
        "%f (mm)  Bmag=%8.4g KG   angle= %g\n",
        (track.P() / units::MeV), Curvature(track) * units::mm, curvaturePlus * units::mm, crtstep / units::mm,
        bmag * toKiloGauss, angle);
// Print("\n");
#endif

  ThreeVector Direction(track.Dx(), track.Dy(), track.Dz());
  ThreeVector PositionNew(0., 0., 0.);
  ThreeVector DirectionNew(0., 0., 0.);

  // char method= '0';
  ThreeVector PositionNewCheck(0., 0., 0.);
  ThreeVector DirectionNewCheck(0., 0., 0.);

  if (useRungeKutta || !fieldConfig->IsFieldUniform()) {
    assert(fieldPropagator);
    fieldPropagator->DoStep(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew, DirectionNew);
#ifdef DEBUG_FIELD
// cross check
#ifndef CHECK_VS_BZ
    ConstFieldHelixStepper stepper(BfieldInitial * toKiloGauss);
    stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNewCheck,
                           DirectionNewCheck);
#else
    double Bz = BfieldInitial[2] * toKiloGauss;
    ConstBzFieldHelixStepper stepper_bz(Bz); //
    stepper_bz.DoStep<ThreeVector, double, int>(Position, Direction, track.Charge(), track.P(), crtstep,
                                                PositionNewCheck, DirectionNewCheck);
#endif

    double posShift = (PositionNew - PositionNewCheck).Mag();
    double dirShift = (DirectionNew - DirectionNewCheck).Mag();

    if (posShift > epsilonRK || dirShift > epsilonRK) {
      std::cout << "*** position/direction shift RK vs. HelixConstBz :" << posShift << " / " << dirShift << "\n";
      printf("*** position/direction shift RK vs. HelixConstBz : %f / %f \n", posShift, dirShift);
      if (verboseDiff) {
        printf("%s End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f\n", " FPH::PiV(1)-RK: ", PositionNew.x(),
               PositionNew.y(), PositionNew.z(), DirectionNew.x(), DirectionNew.y(), DirectionNew.z());
        printf("%s End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f\n", " FPH::PiV(1)-Bz: ", PositionNewCheck.x(),
               PositionNewCheck.y(), PositionNewCheck.z(), DirectionNewCheck.x(), DirectionNewCheck.y(),
               DirectionNewCheck.z());
      }
    }
#endif
// method= 'R';
#ifdef STATS_METHODS
    numRK++;
#endif
  } else {
    // Must agree with values in magneticfield/inc/Units.h
    double Bz             = BfieldInitial[2] * toKiloGauss;
    const bool dominantBz = false;
    if (dominantBz) {
      // Constant field in Z-direction
      ConstBzFieldHelixStepper stepper(Bz); //
      stepper.DoStep<ThreeVector, double, int>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew,
                                               DirectionNew);
// method= 'z';
#ifdef STATS_METHODS
      numHelixZ++;
#endif
    } else {
      // geant::
      double BfieldArr[3] = {BfieldInitial.x() * toKiloGauss, BfieldInitial.y() * toKiloGauss,
                             BfieldInitial.z() * toKiloGauss};
      ConstFieldHelixStepper stepper(BfieldArr);
      stepper.DoStep<double>(Position, Direction, track.Charge(), track.P(), crtstep, PositionNew, DirectionNew);
// method= 'v';
#ifdef STATS_METHODS
      numHelixGen++;
#endif
    }
  }

#ifdef PRINT_FIELD
  // Print(" FPH::PiV(1): Start>", " Pos= %8.5f %8.5f %8.5f  Mom= %8.5f %8.5f %8.5f", Position.x(), Position.y(),
  // Position.z(), Direction.x(), Direction.y(), Direction.z() );
  // Print(" FPH::PiV(1): End>  ", " Pos= %8.5f %8.5f %8.5f  Mom= %8.5f %8.5f %8.5f", PositionNew.x(), PositionNew.y(),
  // PositionNew.z(), DirectionNew.x(), DirectionNew.y(), DirectionNew.z() );

  // printf(" FPH::PiV(1): ");
  printf(" FPH::PiV(1):: ev= %3d trk= %3d %3d %c ", track.Event(), track.Particle(), track.GetNsteps(), method);
  printf("Start> Pos= %8.5f %8.5f %8.5f  Mom= %8.5f %8.5f %8.5f ", Position.x(), Position.y(), Position.z(),
         Direction.x(), Direction.y(), Direction.z());
  printf(" s= %10.6f ang= %7.5f ", crtstep / units::mm, angle);
  printf( // " FPH::PiV(1): "
      "End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f\n", PositionNew.x(), PositionNew.y(), PositionNew.z(),
      DirectionNew.x(), DirectionNew.y(), DirectionNew.z());
#endif

#ifdef STATS_METHODS
  unsigned long modbase = 100;
  if (numTot % modbase < 1) {
    PrintStats();
    if (numTot > 10 * modbase) modbase = 10 * modbase;
  }
#endif

  //  may normalize direction here  // vecCore::math::Normalize(dirnew);
  ThreeVector DirectionUnit = DirectionNew.Unit();
  double posShiftSq         = (PositionNew - Position).Mag2();

  track.SetPosition(PositionNew);
  track.SetDirection(DirectionUnit);
  track.NormalizeFast();

  // Reset relevant variables
  track.SetStatus(kInFlight);
  track.IncreaseStep(crtstep);

  track.DecreasePstep(crtstep);
  if (track.GetPstep() < 1.E-10) {
    track.SetPstep(0);
    track.SetStatus(kPhysics);
  }
  track.DecreaseSnext(crtstep);
  if (track.GetSnext() < 1.E-10) {
    track.SetSnext(0);
  }

  double preSafety = track.GetSafety();
  if (posShiftSq > preSafety * preSafety) {
    track.SetSafety(0);
  } else {
    double posShift = std::sqrt(posShiftSq);
    track.DecreaseSafety(posShift);
    if (track.GetSafety() < 1.E-10) track.SetSafety(0);
  }

#ifdef REPORT_AND_CHECK
  CheckTrack(track, "End of Propagate-In-Volume", 1.0e-5);
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::PropagateInVolume(TrackVec_t &tracks, const double *stepSize, TaskData *td)
{
  // The Vectorized Implementation for Magnetic Field Propagation

  using ThreeVector            = vecgeom::Vector3D<double>;
  constexpr double toKiloGauss = 1.0 / units::kilogauss; // Converts to kilogauss
  bool useRungeKutta           = td->fPropagator->fConfig->fUseRungeKutta;
  const int nTracks            = tracks.size();

  auto fieldConfig = FieldLookup::GetFieldConfig();
  assert(fieldConfig != nullptr);

#if 1 // VECTOR_FIELD_PROPAGATION
  using vecgeom::SOA3D;
  using vecgeom::Vector3D;
  constexpr int Npm = 6;

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

  for (int itr = 0; itr < nTracks; ++itr) {
    Track *pTrack = tracks[itr];

    fltCharge[itr]   = pTrack->Charge();
    momentumMag[itr] = pTrack->P();
    steps[itr]       = stepSize[itr];

    position3D.push_back(pTrack->X(), pTrack->Y(), pTrack->Z());
    direction3D.push_back(pTrack->Dx(), pTrack->Dy(), pTrack->Dz());
  }

  if (!useRungeKutta && fieldConfig->IsFieldUniform()) {
    vecgeom::Vector3D<double> BfieldUniform = fieldConfig->GetUniformFieldValue();
    ConstFieldHelixStepper stepper(BfieldUniform * toKiloGauss);
    // stepper.DoStep<ThreeVector,double,int>(Position,    Direction,  track.Charge(), track.P(), stepSize,
    //                                        PositionNew, DirectionNew);

    // std::cout << "Before Helix stepper - Position addresses: x= " << PositionOut.x() << " y= " << PositionOut.y()
    //           << " z=" << PositionOut.z() << std::endl;

    stepper.DoStepArr<Double_v>(position3D.x(), position3D.y(), position3D.z(), direction3D.x(), direction3D.y(),
                                direction3D.z(), fltCharge, momentumMag, stepSize, PositionOut.x(), PositionOut.y(),
                                PositionOut.z(), DirectionOut.x(), DirectionOut.y(), DirectionOut.z(), nTracks);

    for (int itr = 0; itr < nTracks; ++itr) {
      Track &track = *tracks[itr];

      Vector3D<double> startPosition  = {track.X(), track.Y(), track.Z()};
      Vector3D<double> startDirection = {track.Dx(), track.Dy(), track.Dz()};

#ifdef DEBUG_FIELD
      // Quick Crosscheck against helix stepper
      ThreeVector PositionNew(0., 0., 0.), DirectionNew(0., 0., 0.);

      stepper.DoStep<double>(startPosition, startDirection, track.Charge(), track.P(), stepSize[itr], PositionNew,
                             DirectionNew);

      double posDiff = (PositionNew - PositionOut[itr]).Mag();
      double dirDiff = (DirectionNew - DirectionOut[itr]).Mag();
      if (posDiff > 1.e-6 || dirDiff > 1.e-6) {
        std::cout << "*** position/direction shift HelixStepper scalar vs. vector :" << posDiff << " / " << dirDiff
                  << "\n";
      }
#endif
      Vector3D<double> positionMove = startPosition - PositionOut[itr];

      // Store revised positions and location in original tracks
      track.SetPosition(PositionOut.x(itr), PositionOut.y(itr), PositionOut.z(itr));
      track.SetDirection(DirectionOut.x(itr), DirectionOut.y(itr), DirectionOut.z(itr));
      track.NormalizeFast();

      // Update status, step and safety
      track.SetStatus(kInFlight);
      track.IncreaseStep(stepSize[itr]);

      track.DecreasePstep(stepSize[itr]);
      if (track.GetPstep() < 1.E-10) {
        track.SetPstep(0);
        track.SetStatus(kPhysics);
      }

      track.DecreaseSnext(stepSize[itr]);
      if (track.GetSnext() < 1.E-10) {
        track.SetSnext(0);
      }

      double posShiftSq = positionMove.Mag2();
      double preSafety  = track.GetSafety();
      if (posShiftSq > preSafety * preSafety) {
        track.SetSafety(0);
      } else {
        double posShift = std::sqrt(posShiftSq);
        track.DecreaseSafety(posShift); //  Was crtstep;
        if (track.GetSafety() < 1.0e-10) track.SetSafety(0);
      }
    }
  } else {
    // Prepare for Runge Kutta stepping
    // Choice 3.   Array of FieldTrack
    FieldTrack fldTracksIn[nTracks], fldTracksOut[nTracks];
    for (int itr = 0; itr < nTracks; ++itr) {
      Track *pTrack         = tracks[itr];
      double trackVals[Npm] = {pTrack->X(), pTrack->Y(), pTrack->Z(), pTrack->Px(), pTrack->Py(), pTrack->Pz()};
      fldTracksIn[itr].LoadFromArray(trackVals, Npm);
    }

    auto fieldPropagator = td->fFieldPropagator;
    assert(fieldPropagator);
    auto vectorDriver = fieldPropagator ? fieldPropagator->GetFlexibleIntegrationDriver() : nullptr;
    assert(vectorDriver);

    if (vectorDriver) {
      // Integrate using Runge Kutta method
      vectorDriver->AccurateAdvance(fldTracksIn, steps, fltCharge, fEpsTol, fldTracksOut, nTracks, succeeded);

#ifdef CHECK_VS_SCALAR
      bool checkVsScalar = true;
#ifdef CHECK_VS_HELIX
      const char *diffBanner = "Differences between vector RK vs scalar Helix method:";
#else
      const char *diffBanner = "Differences between vector RK vs scalar RK method:";
#endif
      bool bannerUsed = false;
#endif

      // Store revised positions and location in original tracks
      for (int itr = 0; itr < nTracks; ++itr) {
        Track &track                    = *tracks[itr];
        FieldTrack &fldTrackEnd         = fldTracksOut[itr];
        Vector3D<double> startPosition  = {track.X(), track.Y(), track.Z()};
        Vector3D<double> startDirection = {track.Dx(), track.Dy(), track.Dz()};
        Vector3D<double> endPosition    = {fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2]};

        const double pmag_inv    = 1.0 / track.P();
        ThreeVector endDirVector = pmag_inv * ThreeVector(fldTrackEnd[3], fldTrackEnd[4], fldTrackEnd[5]);

#ifdef CHECK_VS_SCALAR
        double posShift = (startPosition - endPosition).Mag();

        // ---- Perform checks
        ThreeVector endPositionScalar(0., 0., 0.), endDirScalar(0., 0., 0.);
        fieldPropagator->DoStep(startPosition, startDirection, track.Charge(), track.P(), stepSize[itr],
                                endPositionScalar, endDirScalar);

        double posErr = (endPositionScalar - endPosition).Mag();
        double dirErr = (endDirScalar - endDirVector).Mag();
        if (posErr > 1.e-3 * posShift || dirErr > 1.e-6) {
          std::cout << "*** position/direction shift scalar RK vs. vector RK :" << posErr << " / " << dirErr << "\n";
        }

        if (magDiff > perMillion) {
          if (!bannerUsed) {
            std::cerr << diffBanner << std::endl;
            bannerUsed = true;
          }
          std::cerr << "Track " << itr << " has momentum magnitude difference " << magDiff
                    << "  Momentum magnitude @ end = " << std::sqrt(pMag2End) << " vs. start = " << track.P()
                    << std::endl;
          assert(pMag2End > 0.0 && fabs(magDiff) < 0.003 && "ERROR in direction normal.");
        }
#endif

#ifdef DEBUG_FIELD
        // ---- Start verbose print (of selected events/tracks)
        // int maxPartNo = 2, maxEvSlot = 3;
        bool printTrack = false; // = (track.Particle() < maxPartNo) && (track.EventSlot() < maxEvSlot );
        if (printTrack) {
          // Select a few tracks to print ...

          printf(" FPH::PiV(V)/rk: ev= %3d trk= %3d Start> Pos= %8.5f %8.5f %8.5f  Mom= %8.5f %8.5f %8.5f ",
                 track.Event(), track.Particle(), startPosition.x(), startPosition.y(), startPosition.z(),
                 startDirection.x(), startDirection.y(), startDirection.z());
          double angle = std::acos(endDirVector.Dot(startDirection));
          printf(" s= %10.6f ang= %7.5f ", *stepSize / units::mm, angle);
          printf("End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f\n", endPosition.x(), endPosition.y(),
                 endPosition.z(), endDirVector.x(), endDirVector.y(), endDirVector.z());
        }
// ---- End verbose print
#endif

#ifdef CHECK_VS_SCALAR
        // ---- Perform checks
        if (checkVsScalar) {
          bool checkVsHelix = true;
          double curv       = Curvature(track);
          CheckVsScalar(startPosition, startDirection, track.Charge(), track.P(), stepSize[itr], endPosition,
                        endDirVector, curv, itr, td, checkVsHelix);
        }
#endif
        // Update the state of this track
        track.SetPosition(fldTrackEnd[0], fldTrackEnd[1], fldTrackEnd[2]);
        track.SetDirection(endDirVector);
        track.NormalizeFast();
        track.SetStatus(kInFlight);
        track.IncreaseStep(stepSize[itr]);

        track.DecreasePstep(stepSize[itr]);
        if (track.GetPstep() < 1.E-10) {
          track.SetPstep(0);
          track.SetStatus(kPhysics);
        }

        track.DecreaseSnext(stepSize[itr]);
        if (track.GetSnext() < 1.E-10) {
          track.SetSnext(0);
        }

        // Exact update of the safety - using true move (not distance along curve)
        // track.DecreaseSafety(posShift) // More accurate
        double posShiftSq = (startPosition - endPosition).Mag2();
        double preSafety  = track.GetSafety();
        if (posShiftSq > preSafety * preSafety) {
          track.SetSafety(0);
        } else {
          double posShift = std::sqrt(posShiftSq);
          track.DecreaseSafety(posShift); //  Was crtstep;
          if (track.GetSafety() < 1.0e-10) track.SetSafety(0);
        }
      }
    } else {
      // geant::Error( ... );
      std::cerr << "ERROR in FieldPropagationHandler: no Flexible/Vector Integration Driver found." << std::endl;
      exit(1);
    }
  }
#else
  // Placeholder - implemented just as a loop
  for (int itr = 0; itr < nTracks; ++itr)
    PropagateInVolume(*tracks[itr], stepSize[itr], td);
#endif // VECTOR_FIELD_PROPAGATION
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
  // if (track.GetStep() < 1.E-8) td->fNsmall++;
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
void FieldPropagationHandler::PrintStats() const
{
#ifdef STATS_METHODS
  unsigned long nTot = numTot++;
  unsigned long rk = numRK, hZ = numHelixZ, hGen = numHelixGen;
  std::cerr << "Step statistics (field Propagation):  total= " << nTot << " RK = " << rk << "  HelixGen = " << hGen
            << " Helix-Z = " << hZ << std::endl;
#endif
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void FieldPropagationHandler::CheckVsScalar(const vecgeom::Vector3D<double> &startPosition,
                                            const vecgeom::Vector3D<double> &startDirection, double charge,
                                            double momentum, // starting magnitude
                                            double stepSize, const vecgeom::Vector3D<double> &endPosition,
                                            const vecgeom::Vector3D<double> &endDirection, double curvature, int index,
                                            TaskData *td, bool checkVsHelix)
{
  using ThreeVector            = vecgeom::Vector3D<double>;
  constexpr double toKiloGauss = 1.0 / units::kilogauss;
  // static bool bannerUsed  = false;
  auto fieldPropagator = td->fFieldPropagator;
  bool useRungeKutta   = td->fPropagator->fConfig->fUseRungeKutta;
  bool differs         = false;

  // Check against 'scalar' propagation of the same track
  ThreeVector EndPositionScalar(0., 0., 0.), EndDirScalar(0., 0., 0.);

  vecgeom::Vector3D<double> BfieldInitial;
  double bmag;
  FieldLookup::GetFieldValue(startPosition, BfieldInitial, bmag);
  BfieldInitial *= toKiloGauss;
  bmag *= toKiloGauss;

  checkVsHelix = checkVsHelix || (!useRungeKutta); // Cannot check vs RK if it is not available

  if (checkVsHelix) {
    ConstFieldHelixStepper stepper(BfieldInitial);
    stepper.DoStep<double>(startPosition, startDirection, charge, momentum, stepSize, EndPositionScalar, EndDirScalar);
  } else {
    fieldPropagator->DoStep(startPosition, startDirection, charge, momentum, stepSize, EndPositionScalar, EndDirScalar);
  }

  //      checking direction
  ThreeVector diffDir     = endDirection - EndDirScalar;
  double diffDirMag       = diffDir.Mag();
  const double maxDiffMom = 1.0e-4; // 1.5 * fEpsTol; // 10.0 * units::perMillion;
  if (diffDirMag > maxDiffMom) {
    differs = true;
    // if (!bannerUsed) { std::cerr << diffBanner << std::endl; bannerUsed = true; }
    std::ostringstream strDiff;
    strDiff // std::cerr
        << "Track [" << index << "] : direction differs "
        << " by " << diffDir << "  ( mag = " << diffDirMag
        << " ) "
        // << " Direction vector = " << endDirection << "  scalar = " << EndDirScalar
        << " stepSize = " << stepSize << " curv = " << curvature << " B-field = " << BfieldInitial
        << " kilo-Gauss "
        // << " End position= " << endPosition
        << std::endl;
    std::cout << strDiff.str();
    std::cerr << strDiff.str();
  } else {
    //   checking against magnitude of direction difference
    ThreeVector changeDirVector  = endDirection - startDirection;
    ThreeVector changeDirScalar  = EndDirScalar - startDirection;
    double magChangeDirScalar    = changeDirScalar.Mag();
    const double relativeDiffMax = 1.0e-3;

    if (diffDirMag > relativeDiffMax * magChangeDirScalar && (diffDirMag > 1.e-9) && (magChangeDirScalar > 1.e-10)) {
      differs = true;
      // if (!bannerUsed) { std::cerr << diffBanner << std::endl; bannerUsed = true; }
      std::ostringstream strDiff;
      strDiff // std::cerr
          << "Track [" << index << "] : direction CHANGE  has high RELATIVE difference "
          << " by " << diffDir << "  ( mag = " << diffDirMag << " vs |delta Dir scalar| " << magChangeDirScalar
          << "  ) "
          << " DeltaDir/V = " << changeDirVector << "  scalar = " << changeDirScalar << " End position= " << endPosition
          << std::endl;
      std::cout << strDiff.str();
      std::cerr << strDiff.str();
    }
  }

  //      checking position
  ThreeVector diffPos     = endPosition - EndPositionScalar;
  const double maxDiffPos = 1.5 * fEpsTol; // * distanceAlongPath
  if (diffPos.Mag() > maxDiffPos) {
    differs = true;
    // if (!bannerUsed) { std::cerr << diffBanner << std::endl; bannerUsed = true; }
    std::cerr << "Track [" << index << "] : position diff " << diffPos << " mag= " << diffPos.Mag()
              << " "
              // << " Pos: vec = " << endPosition << " scalar = " << EndPositionScalar
              << " move (vector) = " << endPosition - startPosition
              << " its-mag= " << (endPosition - startPosition).Mag() << " stepSize = " << stepSize
              << " move (scalar) = " << EndPositionScalar - startPosition << std::endl;
  }

  if (differs) {
    printf(" FPH::PiV-Start:  Start> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f\n", startPosition.x(),
           startPosition.y(), startPosition.z(), startDirection.x(), startDirection.y(), startDirection.z());
    printf(" FPH::PiV/Vector:   End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f  Delta-p: %6.2g %6.2g %6.2g \n",
           endPosition.x(), endPosition.y(), endPosition.z(), endDirection.x(), endDirection.y(), endDirection.z(),
           endDirection.x() - startDirection.x(), endDirection.y() - startDirection.y(),
           endDirection.z() - startDirection.z());
    printf(" FPH::PiV-1/chk:    End> Pos= %9.6f %9.6f %9.6f  Mom= %9.6f %9.6f %9.6f  Delta-p: %6.2g %6.2g %6.2g \n",
           EndPositionScalar.x(), EndPositionScalar.y(), EndPositionScalar.z(), EndDirScalar.x(), EndDirScalar.y(),
           EndDirScalar.z(), EndDirScalar.x() - startDirection.x(), EndDirScalar.y() - startDirection.y(),
           EndDirScalar.z() - startDirection.z());
  }
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
