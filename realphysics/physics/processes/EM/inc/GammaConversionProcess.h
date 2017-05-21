
#ifndef GAMMACONVERSIONPROCESS_H
#define GAMMACONVERSIONPROCESS_H

#include "EMPhysicsProcess.h"

#include <string>

namespace geantphysics {

class GammaConversionProcess : public EMPhysicsProcess {
public:
  GammaConversionProcess(const std::string &name = "Pair");
};

}        // namespace geantphysics

#endif   // GAMMACONVERSIONPROCESS_H
