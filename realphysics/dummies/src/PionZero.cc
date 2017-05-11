#include "PionZero.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

PionZero* PionZero::Definition() {
  static PionZero instance("pi0", 111, 12, 0.134977, 0);
  return &instance;
}

} // namespace geantphysics
