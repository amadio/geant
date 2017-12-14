
#ifndef GAMMACONVERSIONPROCESS_H
#define GAMMACONVERSIONPROCESS_H

#include "EMPhysicsProcess.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Pre-prepared physics process to describe conversion of photons into e-/e+ pair.
 * @class   GammaConversionProcess
 * @author  M Novak
 * @date    April 2017
 */


class GammaConversionProcess : public EMPhysicsProcess {
public:
  GammaConversionProcess(const std::string &name = "Pair");
};

}        // namespace geantphysics

#endif   // GAMMACONVERSIONPROCESS_H
