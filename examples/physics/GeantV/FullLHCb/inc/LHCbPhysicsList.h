
#ifndef LHCbPHYSICSLIST_H
#define LHCbPHYSICSLIST_H

#include "PhysicsList.h"

#include <string>

namespace lhcbapp {

class LHCbPhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  LHCbPhysicsList(const std::string &name="LHCb-PhysicsList");
  // DTR
  virtual ~LHCbPhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();
};

}           // namespace cmsapp

#endif      // LHCbPHYSICSLIST_H
