
#ifndef CMSPHYSICSLIST_H
#define CMSPHYSICSLIST_H

#include "GeantConfig.h"
#include "Geant/PhysicsList.h"

#include <string>

namespace cmsapp {

class CMSPhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  CMSPhysicsList(const geant::GeantConfig &config, const std::string &name = "CMS-PhysicsList",
                 bool useSamplingTables = false);
  // DTR
  virtual ~CMSPhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

private:
  bool fUseSamplingTables = false;
  bool fVectorized        = false;
  bool fVectorizedMSC     = false;
};

} // namespace cmsapp

#endif // CMSPHYSICSLIST_H
