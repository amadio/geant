
#ifndef LHCbDETECTORCONSTRUCTION_H
#define LHCbDETECTORCONSTRUCTION_H

#include "GeantVDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class RunManager;
  }
}


namespace lhcbapp {

class LHCbDetectorConstruction : public geant::GeantVDetectorConstruction {
public:

  LHCbDetectorConstruction(geant::RunManager *runmgr);

  virtual ~LHCbDetectorConstruction();

  // interface method to define the geometry for the application
  virtual void CreateGeometry();

  void SetGDMLFile(const std::string& gdml) { fGDMLFileName = gdml; }

private:

  std::string  fGDMLFileName;

};

} // namespace lhcbapp

#endif
