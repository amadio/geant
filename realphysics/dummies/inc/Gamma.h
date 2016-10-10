#ifndef GAMMA_H
#define GAMMA_H

#include "Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store gamma static properties.
 * @class   Gamma
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class Gamma : public Particle {
public:
  static Gamma* Definition();

  // copy CTR and assignment operators are deleted
  Gamma(const Gamma&) = delete;
  Gamma& operator=(const Gamma&) = delete;

private:
  Gamma(const std::string &name, int pdgcode, int intcode, double mass, double charge)
  : Particle (name, pdgcode, intcode, mass, charge) {}
};

} // namespace geantphysics

#endif  // GAMMA_H
