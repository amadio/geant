
#include "Geant/PositronAnnihilationProcess.h"

#include "Geant/PositronTo2GammaModel.h"
#include "Geant/Positron.h"
#include "Geant/Gamma.h"

#include "Geant/LightTrack.h"
#include "Geant/PhysicsData.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

PositronAnnihilationProcess::PositronAnnihilationProcess(const std::string &name) : EMPhysicsProcess(name)
{
  // process type is kElectromagnetic in the base EMPhysicsProcess calss
  // set to be a discrete process
  SetIsDiscrete(true);
  // set to be at rest
  SetIsAtRest(true);
  // fill the list of particles that this process can be used to i.e. gamma particle
  AddToListParticlesAlloedToAssigned(Positron::Definition());
  // request to build lambda tables (by default it will be per material)
  //  RequestLambdaTables();
}

void PositronAnnihilationProcess::Initialize()
{
  // model to handle the discrete part pf the interaction (within the EM framework; at-rest part
  // will be handled directly by the process itself in the AtRestDoIt method)
  PositronTo2GammaModel *mod = new PositronTo2GammaModel();
  mod->SetLowEnergyUsageLimit(100. * geant::units::eV);
  mod->SetHighEnergyUsageLimit(100. * geant::units::TeV);
  AddModel(mod);
  // call the EMPhysicsProcess base class init method at the end (after models has been added)
  EMPhysicsProcess::Initialize();
}

double PositronAnnihilationProcess::AverageLifetime(const LightTrack &track) const
{
  if (track.GetKinE() <= 0.) {
    return 0.;
  }
  return GetAVeryLargeValue();
}

int PositronAnnihilationProcess::AtRestDoIt(LightTrack &track, geant::TaskData *td)
{
  int numSecondaries = 2;
  // sample random direction of first photon
  double *rndArray = td->fDblArray;
  td->fRndm->uniform_array(2, rndArray);
  double cost = 2. * rndArray[0] - 1.;
  double sint = std::sqrt((1. - cost) * (1.0 + cost));
  double phi  = geant::units::kTwoPi * rndArray[1];
  // create the 2 secondary particles i.e. the gammas
  numSecondaries = 2;
  // first gamma
  const double xdir       = sint * Math::Cos(phi);
  const double ydir       = sint * Math::Sin(phi);
  LightTrack &gamma1Track = td->fPhysicsData->InsertSecondary();
  gamma1Track.SetDirX(xdir);
  gamma1Track.SetDirY(ydir);
  gamma1Track.SetDirZ(cost);
  gamma1Track.SetKinE(geant::units::kElectronMassC2);
  gamma1Track.SetGVcode(Gamma::Definition()->GetInternalCode()); // gamma GV code
  gamma1Track.SetMass(0.0);
  gamma1Track.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  // second gamma
  LightTrack &gamma2Track = td->fPhysicsData->InsertSecondary();
  gamma2Track.SetDirX(-xdir);
  gamma2Track.SetDirY(-ydir);
  gamma2Track.SetDirZ(-cost);
  gamma2Track.SetKinE(geant::units::kElectronMassC2);
  gamma2Track.SetGVcode(Gamma::Definition()->GetInternalCode()); // gamma GV code
  gamma2Track.SetMass(0.0);
  gamma2Track.SetTrackIndex(track.GetTrackIndex()); // parent Track index
  // kill the primary e+
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);
  //
  return numSecondaries;
}

} // namespace geantphysics
