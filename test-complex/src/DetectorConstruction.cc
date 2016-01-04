
#include "DetectorConstruction.hh"
#include "TGeoManager.h"
#include "TabulatedDataManager.hh"
#include "MaterialConverter.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

TGeoManager * DetectorConstruction::fgGeomMgrRoot= 0 ; // Pointer to the geometry manager   


DetectorConstruction::DetectorConstruction(G4VPhysicalVolume *setWorld) : magField(0){   
      fWorld = setWorld;
/*      char* gdmlFileName = getenv("VP_GEOM_GDML");
      if( gdmlFileName ){
        std::cout << " Creating empty TGeoManager by reading Root geometry from file " << gdmlFileName  << G4endl;
        fgGeomMgrRoot = TGeoManager::Import(gdmlFileName);
      } else {
        std::cout << " Creating empty TGeoManager " << std::endl;
        fgGeomMgrRoot = new TGeoManager();
      }
*/
//      std::cout << " Creating empty TGeoManager " << std::endl;
      fgGeomMgrRoot = new TGeoManager();

      TabulatedDataManager::SetTGeomManager( fgGeomMgrRoot );
      MaterialConverter::SetTGeomManager( fgGeomMgrRoot );
      MaterialConverter::Instance();// call just to initialize

     // initialize magnetic field :: same value as in the prototype
     SetMagField(40.0*CLHEP::kilogauss);
}

DetectorConstruction::~DetectorConstruction(){
  if(magField) delete magField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMagField(G4double fieldValue)
{
  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr
   = G4TransportationManager::GetTransportationManager()->GetFieldManager();

//  if(magField) delete magField;                //delete the existing magn field

  if(fieldValue!=0.)                        // create a new one if non nul
  { magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
    fieldMgr->SetDetectorField(magField);
    fieldMgr->CreateChordFinder(magField);
  } else {
    magField = 0;
    fieldMgr->SetDetectorField(magField);
  }
}

