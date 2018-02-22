#include "KaonPlus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonPlus* KaonPlus::Definition() {
  static KaonPlus instance("K+", 321, 13, 0.493677*geant::units::GeV, geant::units::eplus); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
