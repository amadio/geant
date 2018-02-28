#include "Geant/Gamma.h"

//#include "Geant/PhysicalConstants.h"
//#include "Geant/SystemOfUnits.h"

namespace geantphysics {

Gamma *Gamma::Definition()
{
  static Gamma instance("gamma", 22, 42, 0.0, 0.0);
  return &instance;
}

} // namespace geantphysics
