#include "Geant/Neutron.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

Neutron *Neutron::Definition()
{
  static Neutron instance("neutron", 2112, 3, geant::units::kNeutronMassC2, 0);
  return &instance;
}

} // namespace geantphysics
