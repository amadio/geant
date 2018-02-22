#include "PionPlus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

PionPlus* PionPlus::Definition() {
  static PionPlus instance("pi+", 211, 10, 0.13957*geant::units::GeV, geant::units::eplus); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
