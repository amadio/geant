#ifndef PROTON_H
#define PROTON_H

#include "Particle.h"

namespace geantphysics {
/**
 * @brief   Class(singletone) to store proton static properties.
 * @class   Electron
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class Proton : public Particle {
public:
  static Proton* Definition();

  // copy CTR and assignment operators are deleted
  Proton(const Proton&) = delete;
  Proton& operator=(const Proton&) = delete;

private:
  Proton(const std::string &name, int pdgcode, int intcode, double mass, double charge)
  : Particle (name, pdgcode, intcode, mass, charge) {}
};

}  // namespace geantphysics

#endif  // PROTON_H
