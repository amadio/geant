
#ifndef COMPTONSCATTERINGPROCESS_H
#define COMPTONSCATTERINGPROCESS_H

#include "EMPhysicsProcess.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Pre-prepared physics process to describe incoherent scattering of photons.
 * @class   ComptonScatteringProcess
 * @author  M Novak
 * @date    April 2017
 */


class ComptonScatteringProcess : public EMPhysicsProcess {
public:
  ComptonScatteringProcess(const std::string &name = "Compton");
};

}       // namespace geantphysics

#endif  // COMPTONSCATTERINGPROCESS_H
