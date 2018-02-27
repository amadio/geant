
#ifndef LHCbDETECTORCONSTRUCTION_H
#define LHCbDETECTORCONSTRUCTION_H

#include "Geant/UserDetectorConstruction.h"
#include "Geant/Typedefs.h"
#include "Geant/Config.h"
#include "Geant/Fwd.h"

#include <string>

namespace geant {
  inline namespace GEANT_IMPL_NAMESPACE {
    class RunManager;
  }
}


namespace lhcbapp {

class LHCbDetectorConstruction : public geant::UserDetectorConstruction {
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
