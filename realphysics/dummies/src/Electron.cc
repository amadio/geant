#include "Electron.h"

#include "PhysicalConstants.h"
#include "SystemOfUnits.h"

namespace geantphysics {

Electron* Electron::Definition() {
  static Electron instance("e-",    11, 22, geant::kElectronMassC2, -1.0*geant::eplus);
  return &instance;
}

} // namespace geantphysics
