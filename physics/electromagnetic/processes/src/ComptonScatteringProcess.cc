
#include "ComptonScatteringProcess.h"

#include "Gamma.h"

namespace geantphysics {

ComptonScatteringProcess::ComptonScatteringProcess(const std::string &name) : EMPhysicsProcess(name) {
  // process type is kElectromagnetic in the base EMPhysicsProcess calss
  // set to be a discrete process
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to i.e. gamma particle
  AddToListParticlesAlloedToAssigned(Gamma::Definition());
  // request to build lambda table per-material(by default)
  RequestLambdaTables();
}

}  // namespace geantphysics
