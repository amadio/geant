#include "Geant/Positron.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

Positron* Positron::Definition() {
  static Positron instance("e+",   -11, 23, geant::units::kElectronMassC2,  1.0*geant::units::eplus);
  return &instance;
}

} // namespace geantphysics
