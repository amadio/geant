
#include "CMSDetectorConstruction.h"

#include "GeantRunManager.h"


namespace cmsapp {

CMSDetectorConstruction::CMSDetectorConstruction(geant::GeantRunManager *runmgr)
: geant::GeantVDetectorConstruction(runmgr), fGDMLFileName("cms.gdml") {}


CMSDetectorConstruction::~CMSDetectorConstruction() {}


void CMSDetectorConstruction::CreateGeometry() {
  std::cout<< "  **** LOADING GEOMETRY FROM GDML = " << fGDMLFileName << std::endl;
  LoadGeometry(fGDMLFileName.c_str());
}

} // namespace cmsapp
