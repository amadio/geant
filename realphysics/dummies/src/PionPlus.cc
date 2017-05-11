#include "PionPlus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

PionPlus* PionPlus::Definition() {
  static PionPlus instance("pi+", 211, 10, 0.13957, geant::eplus);
  return &instance;
}

} // namespace geantphysics
