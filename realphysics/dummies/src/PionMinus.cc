#include "PionMinus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

PionMinus* PionMinus::Definition() {
  static PionMinus instance("pi-", -211, 11, 0.13957, -1.0*geant::eplus);
  return &instance;
}

} // namespace geantphysics
