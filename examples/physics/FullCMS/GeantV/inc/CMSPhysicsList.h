
#ifndef CMSPHYSICSLIST_H
#define CMSPHYSICSLIST_H

#include "Geant/PhysicsList.h"

#include <string>

namespace cmsapp {

class CMSPhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  CMSPhysicsList(const std::string &name = "CMS-PhysicsList", bool withAlias = false);
  // DTR
  virtual ~CMSPhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

private:
    bool fWithAlias;
};

} // namespace cmsapp

#endif // CMSPHYSICSLIST_H
