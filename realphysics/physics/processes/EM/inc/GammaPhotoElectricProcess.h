#ifndef GAMMAPHOTOELECTRICPROCESS_H
#define GAMMAPHOTOELECTRICPROCESS_H

#include "EMPhysicsProcess.h"

#include "Gamma.h"

#include <string>

namespace geantphysics {

class GammaPhotoElectricProcess : public EMPhysicsProcess {
public:
  GammaPhotoElectricProcess(const std::string &name = "gPhotoElectric");
};

} // namespace geantphysics

#endif
