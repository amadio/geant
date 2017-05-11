#include "KaonShort.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonShort* KaonShort::Definition() {
  static KaonShort instance("K_S0", 310, 16, 0.497614, 0);
  return &instance;
}

} // namespace geantphysics
