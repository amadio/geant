//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "cmsExpMagneticField.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction(G4int nphi=4, G4int nz=3, G4double density=8.28*g/cm3)
  : _nPhi(nphi), _nZ(nz), _density(density),
    _defaultMaterial(0), _crystalMaterial(0),
    _solidWorld(0),_logicWorld(0),_physWorld(0),
    _solidVols(0),_logicVols(0),_physVols(0),
    _magField(0),_detectorMessenger(0)
{
  //  ComputeCalorParameters();

  // materials
  DefineMaterials();
 
  // create commands for interactive definition of the calorimeter
  _detectorMessenger = new DetectorMessenger(this);

  // magnetic field
  const char* fieldMapFile = getenv("GP_BFIELD_MAP");
  fieldMapFile = (fieldMapFile) ? fieldMapFile : "cmsExp.mag.3_8T";
  _magField = new cmsExpMagneticField();
  _magField->ReadFieldMap( fieldMapFile );
  _magField->SetFieldType("volumebase");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction() {
  delete _detectorMessenger;
  delete _magField;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  //This function illustrates the possible ways to define materials
 
  G4String symbol;             //a=mass of a mole;
  G4double a, z, density;      //z=mean number of protons;  

  //
  // define Elements
  //
  G4Element* elPb = new G4Element("Lead"    ,symbol="Pb", z=82., 207.2*g/mole );
  G4Element* elW  = new G4Element("Tungsten",symbol="W" , z=74., 183.84*g/mole );
  G4Element* elO  = new G4Element("Oxygen"  ,symbol="O" , z=8. ,  15.9994*g/mole );

  //
  // define simple materials from its atomic composition
  //

  // G4int natoms;
  //_crystalMaterial = // new G4Material("FakeStuff",90,120,ecal_density); // _crystalMaterial,ecal_density,90,120);
  _crystalMaterial =  new G4Material("PbWO4", 8.28*g/cm3, 3);
  // _crystalMaterial ->AddElement( elPb, natoms=1 );
  // _crystalMaterial ->AddElement( elW,  natoms=1 );
  // _crystalMaterial ->AddElement( elO,  natoms=4 );
  _crystalMaterial ->AddElement( elPb, 0.4553445);
  _crystalMaterial ->AddElement( elW,  0.4040084);
  _crystalMaterial ->AddElement( elO,  0.1406470);

  //
  // examples of vacuum
  //

  G4Material* Vacuum =
    new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                           kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // print table
  //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  _defaultMaterial = Vacuum;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructCalorimeter()
{

  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // complete the Calor parameters definition
  //ComputeCalorParameters();
   
  //     
  // World
  //
  _solidWorld = new G4Box("World", 9000, 9000, 16000 );
                         
  _logicWorld = new G4LogicalVolume(_solidWorld,           //its solid
                                    _defaultMaterial,      //its material
                                    "WorldLog");           //its name
                                   
  _physWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(),        //at (0,0,0)
                                 _logicWorld,            //its logical volume                                 
                                 "WorldPhys",            //its name
                                 0,                      //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number
  
  //                               
  // Calorimeter
  //  
  const int crystal_n     = _nPhi*_nZ;
  const G4double ecal_zmin  = -3000.;
  const G4double ecal_zmax  =  3000.;

  const G4double ecal_rmin  =  500.;
  const G4double ecal_rmax  =  5000.;
  const G4double ecal_dz    =  0.5*(ecal_zmax-ecal_zmin)/_nZ;
  const G4double ecal_sphi  =     0.;
  const G4double ecal_dphi  =  2.0*M_PI/_nPhi;

  _solidVols = new G4Tubs*[crystal_n];
  _physVols = new G4VPhysicalVolume*[crystal_n]; // reserveNThings<GPVPhysicalVolume>(crystal_n);
  _logicVols = new G4LogicalVolume*[crystal_n]; // reserveNThings<GPLogicalVolume>(crystal_n);

  G4RotationMatrix *idRot = new G4RotationMatrix();/* (1,0,0, */
                                                   /* 0,1,0, */
                                                   /* 0,0,1); */
  
  G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);

  for ( int j = 0; j < _nZ ; ++j ) {
    for ( int i = 0; i < _nPhi ; ++i ) {
      size_t iptr = i+j*_nPhi;
      // G4VSolid *solid = NULL;

      G4ThreeVector ecal_trans(0,0,ecal_zmin+(2.0*j+1.0)*ecal_dz);

      G4cout<<"DetConst.cc: i,j,iptr="<< i <<' '<< j <<' '<< iptr <<' '<< G4endl;

      _solidVols[iptr] = new G4Tubs("tubs", ecal_rmin, ecal_rmax, ecal_dz, ecal_sphi+i*ecal_dphi,ecal_dphi);

      _logicVols[iptr] = new G4LogicalVolume(_solidVols[iptr], _crystalMaterial, "crystal");

      _physVols[iptr] = new G4PVPlacement(idRot, ecal_trans, "crystal", _logicVols[iptr], 
					 _physWorld, //0, 
					 false, 0);
	
      //Set mother
      _physVols[iptr]->SetMotherLogical(_logicWorld);

      //addLogicalVolumePointers( ecal_log[iptr]);
      //addPhysicalVolumePointers( ecal_phy[iptr]);

      _logicVols[iptr]->SetVisAttributes(simpleBoxVisAtt);
    }
  }

  //PrintCalorParameters();     
  
  //                                        
  // Visualization attributes
  //
  _logicWorld->SetVisAttributes (G4VisAttributes::Invisible);

 /*
  // Below are vis attributes that permits someone to test / play 
  // with the interactive expansion / contraction geometry system of the
  // vis/OpenInventor driver :
 {G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  simpleBoxVisAtt->SetVisibility(true);
  delete logicCalor->GetVisAttributes();
  logicCalor->SetVisAttributes(simpleBoxVisAtt);}

 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  logicLayer->SetVisAttributes(atb);}
  
 {G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  atb->SetForceSolid(true);
  logicAbsorber->SetVisAttributes(atb);}
  
 {//Set opacity = 0.2 then transparency = 1 - 0.2 = 0.8
  G4VisAttributes* atb= new G4VisAttributes(G4Colour(0.0,0.0,1.0,0.2));
  atb->SetForceSolid(true);
  logicGap->SetVisAttributes(atb);}
  */

  //
  //always return the physical World
  //
  return _physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4RunManager.hh"

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(ConstructCalorimeter());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
