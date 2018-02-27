
#ifndef ELECTRONIONIZATIONPROCESS_H
#define ELECTRONIONIZATIONPROCESS_H

#include "Geant/EMPhysicsProcess.h"

#include "Geant/Electron.h"
#include "Geant/Positron.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Pre-prepared physics process to describe ionization process of e-/e+.
 * @class   ElectronIonizationProcess
 * @author  M Novak
 * @date    September 2016
 */


class ElectronIonizationProcess : public EMPhysicsProcess {
public:
  ElectronIonizationProcess(const std::string &name = "eIoni");
};

} // namespace geantphysics

#endif
