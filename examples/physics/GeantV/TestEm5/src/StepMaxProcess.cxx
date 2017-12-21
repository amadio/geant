
#include "StepMaxProcess.h"

#include "Electron.h"
#include "Positron.h"
#include "Gamma.h"

//#include "SystemOfUnits.h"

namespace userapplication {

StepMaxProcess::StepMaxProcess(const std::string &name) : geantphysics::PhysicsProcess(name) {
  // set to be a discrete process
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to i.e. gamma particle
  AddToListParticlesAlloedToAssigned(geantphysics::Electron::Definition());
  AddToListParticlesAlloedToAssigned(geantphysics::Positron::Definition());
  AddToListParticlesAlloedToAssigned(geantphysics::Gamma::Definition());
  fMaxStep = geantphysics::PhysicsProcess::GetAVeryLargeValue();
}


StepMaxProcess::~StepMaxProcess() {}

double StepMaxProcess::PostStepLimitationLength(Geant::GeantTrack * /*track*/, Geant::GeantTaskData * /*td*/, bool /*haseloss*/) {
  return fMaxStep;
}


}  // namespace userapplication
