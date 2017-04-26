#ifndef GUGAMMACOMPTONPROCESS_H
#define GUGAMMACOMPTONPROCESS_H

#include "EMPhysicsProcess.h"
#include <string>

namespace geantphysics {

class GUGammaComptonProcess : public EMPhysicsProcess {
public:
  GUGammaComptonProcess(const std::string &name = "gCompton");
};

} // namespace geantphysics

#endif
