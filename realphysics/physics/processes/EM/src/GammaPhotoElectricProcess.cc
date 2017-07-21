#include "GammaPhotoElectricProcess.h"

namespace geantphysics {

GammaPhotoElectricProcess::GammaPhotoElectricProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be an energy loss process (note: loss tables will be built automatically)
  SetType(ProcessType::kElectromagnetic);
  // set to be a continuous-discrete process
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to
  AddToListParticlesAlloedToAssigned(Gamma::Definition());
}

} // namespace geantphysics
