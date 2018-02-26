
#include "PositronAnnihilationProcess.h"

#include "PositronTo2GammaModel.h"
#include "Positron.h"
#include "Gamma.h"

#include "LightTrack.h"
#include "PhysicsData.h"


#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {


PositronAnnihilationProcess::PositronAnnihilationProcess(const std::string &name) : EMPhysicsProcess(name) {
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


void   PositronAnnihilationProcess::Initialize() {
  // model to handle the discrete part pf the interaction (within the EM framework; at-rest part
  // will be handled directly by the process itself in the AtRestDoIt method)
  PositronTo2GammaModel *mod = new PositronTo2GammaModel();
  mod->SetLowEnergyUsageLimit( 100.*geant::units::eV);
  mod->SetHighEnergyUsageLimit(100.*geant::units::TeV);
  AddModel(mod);
  // call the EMPhysicsProcess base class init method at the end (after models has been added)
  EMPhysicsProcess::Initialize();
}

double PositronAnnihilationProcess::AverageLifetime(const LightTrack &track) const {
  if (track.GetKinE()<=0.) {
    return 0.;
  }
  return GetAVeryLargeValue();
}


int    PositronAnnihilationProcess::AtRestDoIt(LightTrack &track, geant::TaskData *td) {
  int numSecondaries = 2;
  // sample random direction of first photon
  double *rndArray = td->fDblArray;
  td->fRndm->uniform_array(2, rndArray);
  double cost = 2.*rndArray[0]-1.;
  double sint = std::sqrt((1.-cost)*(1.0+cost));
  double phi  = geant::units::kTwoPi*rndArray[1];
  // create the 2 secondary partciles i.e. the gammas
  numSecondaries = 2;
  // current capacity of secondary track container
  int curSizeOfSecList = td->fPhysicsData->GetSizeListOfSecondaries();
  // currently used secondary tracks in the container
  int curNumUsedSecs   = td->fPhysicsData->GetNumUsedSecondaries();
  if (curSizeOfSecList-curNumUsedSecs<numSecondaries) {
    td->fPhysicsData->SetSizeListOfSecondaries(2*curSizeOfSecList);
  }
  int secIndx     = curNumUsedSecs;
  curNumUsedSecs += numSecondaries;
  td->fPhysicsData->SetNumUsedSecondaries(curNumUsedSecs);
  std::vector<LightTrack>& sectracks = td->fPhysicsData->GetListOfSecondaries();
  // first gamma
  const double xdir = sint*std::cos(phi);
  const double ydir = sint*std::sin(phi);
  sectracks[secIndx].SetDirX(xdir);
  sectracks[secIndx].SetDirY(ydir);
  sectracks[secIndx].SetDirZ(cost);
  sectracks[secIndx].SetKinE(geant::units::kElectronMassC2);
  sectracks[secIndx].SetGVcode(Gamma::Definition()->GetInternalCode());  // gamma GV code
  sectracks[secIndx].SetMass(0.0);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent Track index
  // first gamma
  ++secIndx;
  sectracks[secIndx].SetDirX(-xdir);
  sectracks[secIndx].SetDirY(-ydir);
  sectracks[secIndx].SetDirZ(-cost);
  sectracks[secIndx].SetKinE(geant::units::kElectronMassC2);
  sectracks[secIndx].SetGVcode(Gamma::Definition()->GetInternalCode());  // gamma GV code
  sectracks[secIndx].SetMass(0.0);
  sectracks[secIndx].SetTrackIndex(track.GetTrackIndex()); // parent Track index
  // kill the primary e+
  track.SetKinE(0.0);
  track.SetTrackStatus(LTrackStatus::kKill);
  //
  return numSecondaries;
}


}  // namespace geantphysics
