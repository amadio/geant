
#ifndef CMSPHYSICSLIST1_H
#define CMSPHYSICSLIST1_H

#include "Geant/PhysicsList.h"

#include <string>

namespace cmsapp {

class CMSPhysicsList1 : public geantphysics::PhysicsList {
public:
  // CTR
  CMSPhysicsList1(const std::string &name = "CMS-PhysicsList1");
  // DTR
  virtual ~CMSPhysicsList1();
  // interface method to assigne physics-process to particles
  virtual void Initialize();
};

} // namespace cmsapp

#endif // CMSPHYSICSLIST1_H
