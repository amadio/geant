#ifndef KAONZERO_H
#define KAONZERO_H

#include "Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store K0 static properties.
 * @class   KaonZero
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class KaonZero : public Particle {
public:
  static KaonZero*  Definition();

  // copy CTR and assignment operators are deleted
  KaonZero(const KaonZero&) = delete;
  KaonZero& operator=(const KaonZero&) = delete;

private:
  KaonZero(const std::string &name, int pdgcode, int intcode, double mass, double charge)
  : Particle (name, pdgcode, intcode, mass, charge) {}
};

} // namespace geantphysics

#endif  // KAONZERO_H
