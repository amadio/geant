
#ifndef LHCbDETECTORCONSTRUCTION_H
#define LHCbDETECTORCONSTRUCTION_H

#include "GeantVDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "GeantFwd.h"

#include <string>

namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class GeantRunManager;
  }
}


namespace lhcbapp {

class LHCbDetectorConstruction : public Geant::GeantVDetectorConstruction {
public:

  LHCbDetectorConstruction(Geant::GeantRunManager *runmgr);

  virtual ~LHCbDetectorConstruction();

  // interface method to define the geometry for the application
  virtual void CreateGeometry();

  void SetGDMLFile(const std::string& gdml) { fGDMLFileName = gdml; }

private:

  std::string  fGDMLFileName;

};

} // namespace lhcbapp

#endif
