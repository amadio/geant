
#include "TestEm5DetectorConstruction.h"

#include "Geant/RunManager.h"

#include "Geant/SystemOfUnits.h"
#include "Geant/Material.h"
#include "Geant/Region.h"

// vecgeom includes
#include "volumes/LogicalVolume.h"
#include "volumes/Box.h"

#include "Geant/Region.h"

namespace userapplication {

TestEm5DetectorConstruction::TestEm5DetectorConstruction(geant::RunManager *runmgr)
    : geant::UserDetectorConstruction(runmgr)
{
  fTargetMatName = "NIST_MAT_Au"; // default target material is NIST Au
  fTargetLogicalVolumeID = -1;
  fTargetRegionIndx      = -1;
  fTargetX               =  19.296*geant::units::um; // set to thickness by default and the half will be taken in ComputeSetup
  fTargetYZ              =  10.*fTargetX;
  fWorldYZ               = 1.2 * fTargetYZ;
  fWorldX                = 1.2 * fTargetX; 
  //
  fTargetMaterial = geantphysics::Material::NISTMaterial(fTargetMatName);
  fWorldMaterial  = geantphysics::Material::NISTMaterial("NIST_MAT_Galactic");
  //
  fGammaCut    = 13.00 * geant::units::um;
  fElectronCut = 1.12 * geant::units::um;
  fPositronCut = fElectronCut;
}

TestEm5DetectorConstruction::~TestEm5DetectorConstruction()
{ /* nothing to do */
}

// NOTE: custom materials are created only for demonstration but when simulation results are supposed to be compared to
//       the corresponding Geant4 simulation we strongly suggest to use pre-defined NIST materials on both sides!
// we create only 3 custom materials that can be selected by the user by name; all pre-defined NIST materials are
// available and the user can select them by name or the user can add more custom materials to this method
void TestEm5DetectorConstruction::CreateMaterials()
{
  int z, ncomponents, natoms;
  double a, density, massfraction;
  //
  // create simple (single element) custom material Gold:
  new geantphysics::Material("Gold", z = 79, a = 196.97 * geant::units::g / geant::units::mole,
                             density = 19.32 * geant::units::g / geant::units::cm3);
  //
  // create complex custom material Water (adding elements with number of atoms in the molecule):
  // 1. create the custom H and O elements
  geantphysics::Element *elH =
      new geantphysics::Element("Hydrogen", "H", z = 1, a = 1.01 * geant::units::g / geant::units::mole);
  geantphysics::Element *elO =
      new geantphysics::Element("Oxygen", "O", z = 8, a = 16.00 * geant::units::g / geant::units::mole);
  // 2. create the material object by specifying density and number of components to be added
  geantphysics::Material *matH2O =
      new geantphysics::Material("Water", density = 1.000 * geant::units::g / geant::units::cm3, ncomponents = 2);
  matH2O->AddElement(elH, natoms = 2);
  matH2O->AddElement(elO, natoms = 1);
  // 3. add the elements(components) to the material by specifying number of atoms from the given element
  //
  // create complex custom material (adding components by fractional mass) Air:
  // 1. use NIST element N and the previously created custom element elO
  // 2. create the material object by specifying density and number of components to be added
  geantphysics::Material *matAir =
      new geantphysics::Material("Air", density = 1.290 * geant::units::mg / geant::units::cm3, ncomponents = 2);
  matAir->AddElement(geantphysics::Element::NISTElement(z = 7), massfraction = 0.7);
  matAir->AddElement(elO, massfraction = 0.3);
}

void TestEm5DetectorConstruction::ComputeSetup()
{
  fTargetMaterial = geantphysics::Material::NISTMaterial(fTargetMatName);

  fTargetYZ  =  10.*fTargetX;
  fTargetYZ *=  0.5;
  fTargetX  *=  0.5;
  fWorldYZ   =  1.2 * fTargetYZ;
  fWorldX    =  1.2 * fTargetX;
}

void TestEm5DetectorConstruction::CreateGeometry()
{
  //  LoadGeometry("/Users/mnovak/opt/Data/GV/hanson_19.root");
  //  return;
  // compute geometry setup parameters
  ComputeSetup();
  // create one region that will be assigned to all logical volume (i.e. enough to set for the world)
  vecgeom::Region *aRegion = new vecgeom::Region("UserDefinedRegion", true, fGammaCut, fElectronCut, fPositronCut);
  //
  // create geometry
  vecgeom::UnplacedBox *world  = new vecgeom::UnplacedBox(fWorldX, fWorldYZ, fWorldYZ);
  vecgeom::UnplacedBox *target = new vecgeom::UnplacedBox(fTargetX, fTargetYZ, fTargetYZ);

  // create the corresponding logical volumes
  vecgeom::LogicalVolume *logicWorld  = new vecgeom::LogicalVolume("world", world);
  vecgeom::LogicalVolume *logicTarget = new vecgeom::LogicalVolume("target", target);
  // set region
  logicWorld->SetRegion(aRegion);
  logicTarget->SetRegion(aRegion);
  // set target logical colume ID and target region index (TestEm5 application will require it)
  fTargetLogicalVolumeID = logicTarget->id();
  fTargetRegionIndx      = aRegion->GetIndex();
  // set materials
  logicWorld->SetMaterialPtr((void *)fWorldMaterial);
  logicTarget->SetMaterialPtr((void *)fTargetMaterial);
  // place target into world
  vecgeom::Transformation3D placement(0.0, 0, 0);
  logicWorld->PlaceDaughter("target", logicTarget, &placement);
  // set world and close geometry
  vecgeom::GeoManager::Instance().SetWorld(logicWorld->Place());
  vecgeom::GeoManager::Instance().CloseGeometry();
}

} // namespace userapplication
