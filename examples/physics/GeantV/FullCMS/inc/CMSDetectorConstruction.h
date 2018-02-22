
#ifndef CMSDETECTORCONSTRUCTION_H
#define CMSDETECTORCONSTRUCTION_H

#include "GeantVDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class GeantRunManager;
  }
}


namespace cmsapp {

class CMSDetectorConstruction : public geant::GeantVDetectorConstruction {
public:

  CMSDetectorConstruction(geant::GeantRunManager *runmgr);

  virtual ~CMSDetectorConstruction();

  // interface method to define the geometry for the application
  virtual void CreateGeometry();

  void SetGDMLFile(const std::string& gdml) { fGDMLFileName = gdml; }

private:

  std::string  fGDMLFileName;

};

} // namespace cmsapp

#endif
