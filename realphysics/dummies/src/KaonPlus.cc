#include "KaonPlus.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

KaonPlus* KaonPlus::Definition() {
  static KaonPlus instance("K+", 321, 13, 0.493677, geant::eplus);
  return &instance;
}

} // namespace geantphysics
