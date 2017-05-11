#include "KaonZero.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonZero* KaonZero::Definition() {
  static KaonZero instance("K0", 311, 17, 0.497614, 0);
  return &instance;
}

} // namespace geantphysics
