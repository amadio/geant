
#include "LHCbDetectorConstruction.h"

#include "GeantRunManager.h"


namespace lhcbapp {

LHCbDetectorConstruction::LHCbDetectorConstruction(Geant::GeantRunManager *runmgr)
: Geant::GeantVDetectorConstruction(runmgr), fGDMLFileName("LHCb_201603.gdml") {}


LHCbDetectorConstruction::~LHCbDetectorConstruction() {}


void LHCbDetectorConstruction::CreateGeometry() {
  std::cout<< "  **** LOADING GEOMETRY FROM GDML = " << fGDMLFileName << std::endl;
  LoadGeometry(fGDMLFileName.c_str());
}

} // namespace lhcbapp
