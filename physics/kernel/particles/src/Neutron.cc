#include "Neutron.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

Neutron* Neutron::Definition() {
  static Neutron instance("neutron", 2112, 3, geant::kNeutronMassC2, 0);
  return &instance;
}

} // namespace geantphysics
