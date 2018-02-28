#ifndef ELECTRON_H
#define ELECTRON_H

#include "Geant/Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store electron static properties.
 * @class   Electron
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class Electron : public Particle {
public:
  static Electron *Definition();

  // copy CTR and assignment operators are deleted
  Electron(const Electron &) = delete;
  Electron &operator=(const Electron &) = delete;

private:
  Electron(const std::string &name, int pdgcode, int intcode, double mass, double charge)
      : Particle(name, pdgcode, intcode, mass, charge)
  {
  }
};

} // namespace geantphysics

#endif // ELECTRON_H
