#include "Proton.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

Proton* Proton::Definition() {
  static Proton instance("proton", 2212, 1, geant::kProtonMassC2,  1.0*geant::eplus);
  return &instance;
}

} // namespace geantphysics
