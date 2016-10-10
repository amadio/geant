
#ifndef ELECTRONBREMSSTRAHLUNGPROCESS_H
#define ELECTRONBREMSSTRAHLUNGPROCESS_H

#include "EMPhysicsProcess.h"

#include "Electron.h"
#include "Positron.h"

#include <string>

namespace geantphysics {

class ElectronBremsstrahlungProcess : public EMPhysicsProcess {
public:
  // CTR
  ElectronBremsstrahlungProcess(const std::string &name = "eBrem");
};

} // namespace geantphysics

#endif
