
#include "ElectronBremsstrahlungProcess.h"

namespace geantphysics {

ElectronBremsstrahlungProcess::ElectronBremsstrahlungProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be an energy loss process (note: loss tables will be built automatically)
  SetType(ProcessType::kEnergyLoss);
  // set to be a continuous-discrete process
  SetIsContinuous(true);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to (note: models need to set either to be for e- or e+)
  AddToListParticlesAlloedToAssigned(Electron::Definition());
  AddToListParticlesAlloedToAssigned(Positron::Definition());
  // request to build lambda table per-material-cuts (per-material by default)
  RequestLambdaTables(false);
}

} // namespace geantphysics
