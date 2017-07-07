
#include "ElasticScatteringProcess.h"

namespace geantphysics {

ElasticScatteringProcess::ElasticScatteringProcess(const std::string &name)
: HadronicProcess(name) {
  // set process type to be an elastic scattering process
  SetType(HadronicProcessType::kElastic);
  // set to be a continuous-discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to
  AddToListParticlesAlloedToAssigned(Proton::Definition());
  AddToListParticlesAlloedToAssigned(Neutron::Definition());
}

} // namespace geantphysics
