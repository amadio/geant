#ifndef GAMMAPHOTOELECTRICPROCESS_H
#define GAMMAPHOTOELECTRICPROCESS_H

#include "Geant/EMPhysicsProcess.h"

#include "Geant/Gamma.h"

#include <string>

namespace geantphysics {

class GammaPhotoElectricProcess : public EMPhysicsProcess {
public:
  GammaPhotoElectricProcess(const std::string &name = "gPhotoElectric");
};

} // namespace geantphysics

#endif
