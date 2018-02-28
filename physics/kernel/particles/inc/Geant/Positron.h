#ifndef POSITRON_H
#define POSITRON_H

#include "Geant/Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store positron static properties.
 * @class   Electron
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class Positron : public Particle {
public:
  static Positron *Definition();

  // copy CTR and assignment operators are deleted
  Positron(const Positron &) = delete;
  Positron &operator=(const Positron &) = delete;

private:
  Positron(const std::string &name, int pdgcode, int intcode, double mass, double charge)
      : Particle(name, pdgcode, intcode, mass, charge)
  {
  }
};

} // namespace geantphysics

#endif // POSITRON_H
