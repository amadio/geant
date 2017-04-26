#ifndef GUGAMMACONVERSIONPROCESS_H
#define GUGAMMACONVERSIONPROCESS_H

#include "EMPhysicsProcess.h"
#include <string>

namespace geantphysics {

class GUGammaConversionProcess : public EMPhysicsProcess {
public:
  GUGammaConversionProcess(const std::string &name = "gConversion");
};

} // namespace geantphysics

#endif 
