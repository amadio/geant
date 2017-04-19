
#ifndef COMPTONSCATTERINGPROCESS_H
#define COMPTONSCATTERINGPROCESS_H

#include "EMPhysicsProcess.h"

#include <string>

namespace geantphysics {

class ComptonScatteringProcess : public EMPhysicsProcess {
public:
  ComptonScatteringProcess(const std::string &name = "Compton");
};

}       // namespace geantphysics

#endif  // COMPTONSCATTERINGPROCESS_H
