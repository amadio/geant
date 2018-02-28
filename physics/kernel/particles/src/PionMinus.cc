#include "Geant/PionMinus.h"

#include "Geant/PhysicalConstants.h"
#include "Geant/SystemOfUnits.h"

namespace geantphysics {

PionMinus *PionMinus::Definition()
{
  static PionMinus instance("pi-", -211, 11, 0.13957 * geant::units::GeV,
                            -1.0 * geant::units::eplus); // mass value taken from Geant4 10.3
  return &instance;
}

} // namespace geantphysics
