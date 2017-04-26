#ifndef GUGAMMAPHOTONELECTRICPROCESS_H
#define GUGAMMAPHOTONELECTRICPROCESS_H

#include "EMPhysicsProcess.h"
#include <string>

namespace geantphysics {

class GUGammaPhotoElectricProcess : public EMPhysicsProcess {
public:
  GUGammaPhotoElectricProcess(const std::string &name = "gPhotoElectric");
};

} // namespace geantphysics

#endif
