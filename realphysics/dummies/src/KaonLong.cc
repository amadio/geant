#include "KaonLong.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonLong* KaonLong::Definition() {
  static KaonLong instance("K_L0", 130, 15, 0.497614*geant::GeV, 0); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
