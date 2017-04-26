#include "GUGammaConversionProcess.h"
#include "Gamma.h"

namespace geantphysics {

GUGammaConversionProcess::GUGammaConversionProcess(const std::string &name)
  : EMPhysicsProcess(name) 
{
  // set process type to be an electromagnetic process
  SetType(ProcessType::kElectromagnetic);
  // set to be a discrete process
  SetIsContinuous(false);
  SetIsDiscrete(true);
  // fill the list of particles that this process can be used  
  AddToListParticlesAlloedToAssigned(Gamma::Definition());
}

} // namespace geantphysics
