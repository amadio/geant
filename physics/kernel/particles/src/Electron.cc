#include "Geant/Electron.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

Electron *Electron::Definition()
{
  static Electron instance("e-", 11, 22, geant::units::kElectronMassC2, -1.0 * geant::units::eplus);
  return &instance;
}

} // namespace geantphysics
