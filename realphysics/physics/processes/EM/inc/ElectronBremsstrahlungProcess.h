
#ifndef ELECTRONBREMSSTRAHLUNGPROCESS_H
#define ELECTRONBREMSSTRAHLUNGPROCESS_H

#include "EMPhysicsProcess.h"

#include "Electron.h"
#include "Positron.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Pre-prepared physics process to describe bremsstrahlung photon emission of e-/e+.
 * @class   ElectronBremsstrahlungProcess
 * @author  M Novak
 * @date    September 2016
 */


class ElectronBremsstrahlungProcess : public EMPhysicsProcess {
public:
  // CTR
  ElectronBremsstrahlungProcess(const std::string &name = "eBrem");
};

} // namespace geantphysics

#endif
