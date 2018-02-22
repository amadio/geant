
#ifndef CMSDETECTORCONSTRUCTION_H
#define CMSDETECTORCONSTRUCTION_H

#include "UserDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class RunManager;
  }
}


namespace cmsapp {

class CMSDetectorConstruction : public geant::UserDetectorConstruction {
public:

  CMSDetectorConstruction(geant::RunManager *runmgr);

  virtual ~CMSDetectorConstruction();

  // interface method to define the geometry for the application
  virtual void CreateGeometry();

  void SetGDMLFile(const std::string& gdml) { fGDMLFileName = gdml; }

private:

  std::string  fGDMLFileName;

};

} // namespace cmsapp

#endif
