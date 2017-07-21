#include "GammaPhotoElectricProcess.h"

namespace geantphysics {

GammaPhotoElectricProcess::GammaPhotoElectricProcess(const std::string &name)
: EMPhysicsProcess(name) {
  // set process type to be an electromagnetic process
  SetType(ProcessType::kElectromagnetic);
  // set to be a discrete process
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used to
  AddToListParticlesAlloedToAssigned(Gamma::Definition());
}

} // namespace geantphysics
