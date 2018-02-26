#include "KaonShort.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

KaonShort* KaonShort::Definition() {
  static KaonShort instance("K_S0", 310, 16, 0.497614*geant::units::GeV, 0); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
