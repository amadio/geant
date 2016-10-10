
#ifndef ELECTRONIONIZATIONPROCESS_H
#define ELECTRONIONIZATIONPROCESS_H

#include "EMPhysicsProcess.h"

#include "Electron.h"
#include "Positron.h"

#include <string>

namespace geantphysics {

class ElectronIonizationProcess : public EMPhysicsProcess {
public:
  ElectronIonizationProcess(const std::string &name = "eIoni");
};

} // namespace geantphysics

#endif
