
#include "TestEm3DetectorConstruction.h"

#include "GeantVDetectorConstruction.h"
#include "GeantRunManager.h"

#include <iostream>
#include <vector>

// Material includes
#include "Material.h"
#include "Element.h"
#include "MaterialProperties.h"
//#include "NISTElementData.h"
#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

// Region and productionCut includes
#include "Region.h"
#include "PhysicsParameters.h"
#include "MaterialCuts.h"
/////////////////////////////////
// VecGeom includes
#include "management/GeoManager.h"
#include "volumes/Box.h"
#include "volumes/LogicalVolume.h"
/////////////////////////////////

namespace userapplication {

TestEm3DetectorConstruction::TestEm3DetectorConstruction(geant::GeantRunManager *runmgr) : geant::GeantVDetectorConstruction(runmgr) {
  fNumberOfAbsorbers        = 2;
  fAbsorberThicknesses[0]   = 2.3*geant::units::mm;
  fAbsorberThicknesses[1]   = 5.7*geant::units::mm;
  fNumberOfLayers           = 50;
  fCaloSizeYZ               = 40.*geant::units::cm;
  //
  ComputeCalorimeter();
  //
  fWorldMaterialName        = "Galactic";
  fAbsorberMaterialNames[0] = "NIST_MAT_Pb";
  fAbsorberMaterialNames[1] = "NIST_MAT_lAr";
  //
  fProductionCut            = 0.7*geant::units::mm;
  //
  for (int i=0; i<gMaxNumAbsorbers; ++i) {
    fAbsorberLogicVolIDs[i] = -1;
    fAbsorberMaterials[i]   = nullptr;
  }
  fWorldMaterial            = nullptr;
}


TestEm3DetectorConstruction::~TestEm3DetectorConstruction() { /* nothing to do */}



void   TestEm3DetectorConstruction::SetNumberOfAbsorbersPerLayer(int numabs) {
  if (numabs>gMaxNumAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::SetNumberOfAbsorbersPerLayer() \n"
                "    Number of absorbers per layer " << numabs << " > maximum " << gMaxNumAbsorbers
             << std::endl;
    exit(-1);
  }
  fNumberOfAbsorbers = numabs;
}


void   TestEm3DetectorConstruction::SetAbsorberThickness(int absindx, double thick) {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::SetAbsorberThickness() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  fAbsorberThicknesses[absindx] = thick;
}


double TestEm3DetectorConstruction::GetAbsorberThickness(int absindx) const {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::GetAbsorberThickness() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  return fAbsorberThicknesses[absindx];
}


void  TestEm3DetectorConstruction::SetAbsorberMaterialName(int absindx, const std::string &matname) {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::SetAbsorberMaterialName() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  fAbsorberMaterialNames[absindx] = matname;
}


void   TestEm3DetectorConstruction::SetAbsorberMaterial(int absindx, const std::string &matname) {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::SetAbsorberMaterial() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  // will report proper error msg if the material cannot be found or created
  geantphysics::Material *mat = geantphysics::Material::NISTMaterial(matname);
  fAbsorberMaterials[absindx] = mat;
}


const geantphysics::Material* TestEm3DetectorConstruction::GetAbsorberMaterial(int absindx) const {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::GetAbsorberMaterial() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  return fAbsorberMaterials[absindx];
}


int    TestEm3DetectorConstruction::GetAbsorberLogicalVolumeID(int absindx) const {
  if (absindx>=fNumberOfAbsorbers) {
    std::cerr<< "  **** ERROR: TestEm3DetectorConstruction::GetAbsorberLogicalVolumeID() \n"
                "    Unknown absorber index = " << absindx << " ( should be <= " << fNumberOfAbsorbers-1 << " )"
             << std::endl;
    exit(-1);
  }
  return fAbsorberLogicVolIDs[absindx];
}



void   TestEm3DetectorConstruction::ComputeCalorimeter() {
  fLayerThickness = 0.0;
  for (int i=0; i<fNumberOfAbsorbers; ++i) {
    fLayerThickness += fAbsorberThicknesses[i];
  }
  fCaloSizeX   = fNumberOfLayers*fLayerThickness;
  fWorldSizeYZ = 1.2*fCaloSizeYZ;
  fWorldSizeX  = 1.2*fCaloSizeX;
}

void   TestEm3DetectorConstruction::DetectorInfo() {
  std::cout<< "\n ========================    Detector Info   ========================================  "<< std::endl;
  std::cout<< "     Calorimeter is " << fNumberOfLayers << " layers of [ ";
           for (int i=0; i<fNumberOfAbsorbers; ++i) {
             std::cerr<< fAbsorberThicknesses[i]/geant::units::mm << "  mm "
                      << fAbsorberMaterialNames[i];
             if (i<fNumberOfAbsorbers-1) std::cerr<<" + ";
           }
           std::cout<< " ]"<<std::endl;
  std::cout<<" ====================================================================================  \n"<<std::endl;
}


void   TestEm3DetectorConstruction::CreateMaterials() {
  geantphysics::Element* elH   = geantphysics::Element::NISTElement(1);
  geantphysics::Element* elC   = geantphysics::Element::NISTElement(6);
  geantphysics::Element* elPb  = geantphysics::Element::NISTElement(82);
  //
  // create material Galactic
  double density      = geant::units::kUniverseMeanDensity;
  double pressure     = 3.e-18*geant::units::pascal;
  double temperature  = 2.73*geant::units::kelvin;
  double z            = 1.;
  double a            = 1.008*geant::units::g/geant::units::mole;
  new geantphysics::Material("Galactic", z, a, density, geantphysics::MaterialState::kStateGas, temperature, pressure);
  //
  // create material Lead (by mass fraction)
  double massfraction = 1.;
  int ncomponents     = 1;
  density             = 11.35*geant::units::g/geant::units::cm3;
  geantphysics::Material *matLead = new geantphysics::Material("Lead", density, ncomponents);
  matLead->AddElement(elPb, massfraction=1.0);
  //
  // create Scintillator (by atom count)
  density             = 1.032*geant::units::g/geant::units::cm3;
  ncomponents         = 2;
  int natoms          = 1;
  geantphysics::Material *matSci = new geantphysics::Material("Scintillator", density, ncomponents);
  matSci->AddElement(elC, natoms =  9);
  matSci->AddElement(elH, natoms = 10);
  //
  // create liquidArgon (from material NIST_MAT_lAr with different density by mass fraction)
  geantphysics::Material *matlAr = geantphysics::Material::NISTMaterial("NIST_MAT_lAr");
  density             = 1.390*geant::units::g/geant::units::cm3;
  ncomponents         = 1;
  massfraction        = 1.;
  geantphysics::Material *lArEm3 = new geantphysics::Material("liquidArgon", density, ncomponents);
  lArEm3->AddMaterial(matlAr, massfraction);
}


void TestEm3DetectorConstruction::CreateGeometry() {
  // first set the materials that have been defined by their name
  fWorldMaterial = geantphysics::Material::NISTMaterial(fWorldMaterialName);
  for (int i=0; i<fNumberOfAbsorbers; ++i) {
    SetAbsorberMaterial(i,fAbsorberMaterialNames[i]);
  }
  // create one region (with the production cut in length)
  // all logical volume will be assigned to this region (i.e. enough to set for the world)
  vecgeom::Region *aRegion = new vecgeom::Region("Region", true, fProductionCut, fProductionCut, fProductionCut);
  //
  // create geometry
  vecgeom::UnplacedBox *world  = new vecgeom::UnplacedBox(0.5*fWorldSizeX,0.5*fWorldSizeYZ,0.5*fWorldSizeYZ);
  // create the corresponding logical volume
  vecgeom::LogicalVolume *logicWorld  = new vecgeom::LogicalVolume("world",world);
  // set material pointer
  logicWorld->SetMaterialPtr(fWorldMaterial);
  //
  // create layer
  vecgeom::UnplacedBox *layer = new vecgeom::UnplacedBox("layer",0.5*fLayerThickness,0.5*fCaloSizeYZ,0.5*fCaloSizeYZ);
  // create the corresponding logical volume
  vecgeom::LogicalVolume *logicLayer  = new vecgeom::LogicalVolume("layer",layer);
  // set material pointer
  logicLayer->SetMaterialPtr(fWorldMaterial);
  //
  // create all absorbers and place them into the layer
  char name[512];
  double xcenter = 0.;
  double xstart  = -0.5*fLayerThickness;
  for (int i=0; i<fNumberOfAbsorbers; ++i) {
    // create absorber i
    sprintf(name,"abs_%d",i);
    double sizeX  = 0.5*fAbsorberThicknesses[i];
    double sizeYZ = 0.5*fCaloSizeYZ;
    vecgeom::UnplacedBox *absor = new vecgeom::UnplacedBox(name,sizeX,sizeYZ,sizeYZ);
    // create the corresponding logical volume
    vecgeom::LogicalVolume *logicAbsor =  new vecgeom::LogicalVolume(name,absor);
    // set material
    logicAbsor->SetMaterialPtr(fAbsorberMaterials[i]);
    // get the logical volume ID
    fAbsorberLogicVolIDs[i] = logicAbsor->id();
    // compute placement and place the absorber into the layer
    xcenter = xstart+sizeX;
    xstart += fAbsorberThicknesses[i];
    vecgeom::Transformation3D *place = new vecgeom::Transformation3D(xcenter, 0, 0, 0, 0, 0);
    logicLayer->PlaceDaughter(name, logicAbsor, place);
  }
  //
  // create calorimeter and place fNumberOfLayers into the calorimeter
  vecgeom::UnplacedBox *calo = new vecgeom::UnplacedBox("calo",0.5*fCaloSizeX,0.5*fCaloSizeYZ,0.5*fCaloSizeYZ);
  // create the corresponding logical volume
  vecgeom::LogicalVolume *logicCalo  = new vecgeom::LogicalVolume("calo",calo);
  // set material pointer
  logicCalo->SetMaterialPtr(fWorldMaterial);
  xcenter = 0.;
  xstart  = -0.5*fCaloSizeX;
  for (int i=0; i<fNumberOfLayers; ++i) {
    sprintf(name,"layer_%d",i);
    xcenter = xstart+0.5*fLayerThickness;
    xstart += fLayerThickness;
    vecgeom::Transformation3D *place = new vecgeom::Transformation3D(xcenter, 0, 0, 0, 0, 0);
    logicCalo->PlaceDaughter(name, logicLayer, place);
  }
  // place the calorimeter into the world
  vecgeom::Transformation3D *place = new vecgeom::Transformation3D(0, 0, 0, 0, 0, 0);
  logicWorld->PlaceDaughter("calorimeter", logicCalo, place);
  //
  // set the world logical volume, assigned it to the created region and close the geometry
  vecgeom::VPlacedVolume *w = logicWorld->Place();
  vecgeom::GeoManager::Instance().SetWorld(w);
  // set region for the world logical volume (will be set for all daughter volumes)
  logicWorld->SetRegion(aRegion);
  vecgeom::GeoManager::Instance().CloseGeometry();
}


}    // namespace userapplication
