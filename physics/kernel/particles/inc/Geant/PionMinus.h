#ifndef PIONMINUS_H
#define PIONMINUS_H

#include "Geant/Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store pi- static properties.
 * @class   PionMinus
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class PionMinus : public Particle {
public:
  static PionMinus*  Definition();

  // copy CTR and assignment operators are deleted
  PionMinus(const PionMinus&) = delete;
  PionMinus& operator=(const PionMinus&) = delete;

private:
  PionMinus(const std::string &name, int pdgcode, int intcode, double mass, double charge)
  : Particle (name, pdgcode, intcode, mass, charge) {}
};

} // namespace geantphysics

#endif  // PIONMINUS_H
