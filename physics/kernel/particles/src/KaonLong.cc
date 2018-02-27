#include "Geant/KaonLong.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

KaonLong* KaonLong::Definition() {
  static KaonLong instance("K_L0", 130, 15, 0.497614*geant::units::GeV, 0); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
