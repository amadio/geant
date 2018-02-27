
#include "CMSDetectorConstruction.h"

#include "Geant/RunManager.h"


namespace cmsapp {

CMSDetectorConstruction::CMSDetectorConstruction(geant::RunManager *runmgr)
: geant::UserDetectorConstruction(runmgr), fGDMLFileName("cms.gdml") {}


CMSDetectorConstruction::~CMSDetectorConstruction() {}


void CMSDetectorConstruction::CreateGeometry() {
  std::cout<< "  **** LOADING GEOMETRY FROM GDML = " << fGDMLFileName << std::endl;
  LoadGeometry(fGDMLFileName.c_str());
}

} // namespace cmsapp
