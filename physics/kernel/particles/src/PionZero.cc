#include "PionZero.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

PionZero* PionZero::Definition() {
  static PionZero instance("pi0", 111, 12, 0.134977*geant::units::GeV, 0); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
