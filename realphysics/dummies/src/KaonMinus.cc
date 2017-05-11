#include "KaonMinus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonMinus* KaonMinus::Definition() {
  static KaonMinus instance("K-", -321, 14, 0.493677, -1.0*geant::eplus);
  return &instance;
}

} // namespace geantphysics
