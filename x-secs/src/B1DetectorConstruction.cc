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
// $Id$
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1SteppingAction.hh"
#include "B1materials.hh"
   // use of stepping action to set the accounting volume

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

using namespace CLHEP;

G4double MaterialPosition[98][3]; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  G4bool checkOverlaps = true;
  // Envelope parameters
  //
  const G4double kLattice = 15*cm;
  const G4double kRadius = 5*cm;
  const G4int nballs = 10;
  G4double world_sizeXY = kLattice*(nballs-1)+2*(kLattice-kRadius);
  G4double world_sizeZ = 2*(kLattice-kRadius);
  //
  // Create the empty world
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
  //  G4cout << world_mat;
  
  G4Box* solidWorld =    
     new G4Box("World",                       //its name
	       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
  
  G4LogicalVolume* logicWorld =                         
     new G4LogicalVolume(solidWorld,          //its solid
			 world_mat,           //its material
			 "World");            //its name
  
  G4VPhysicalVolume* physWorld = 
     new G4PVPlacement(0,                     //no rotation
		       G4ThreeVector(),       //at (0,0,0)
		       logicWorld,            //its logical volume
		       "World",               //its name
		       0,                     //its mother  volume
		       false,                 //no boolean operation
		       0,                     //copy number
		       checkOverlaps);        //overlaps checking


  // 
  // Now the balls!
  //

  B1SteppingAction* steppingAction = B1SteppingAction::Instance(); 
  for(G4int i=0; i<nmaterials; ++i) {
     G4Material *mat = nist->FindOrBuildMaterial(materialVec[i]);
     //     G4cout << *mat;
     G4Sphere *solidSphere = new G4Sphere(materialVec[i],             // name
					  0.,kRadius,              // rmin rmax
					  0., 2*pi,                // phi range
					  0., pi);                 // phi range
     G4LogicalVolume* logicSphere =                         
	new G4LogicalVolume(solidSphere,            //its solid
			    mat,                    //its material
			    materialVec[i]);           //its name
  
     G4double xplace = kLattice-kRadius+(i%nballs)*kLattice-0.5*world_sizeXY;
     G4double yplace = kLattice-kRadius+(i/nballs)*kLattice-0.5*world_sizeXY;
     //     printf("%f %f\n",xplace,yplace);
     new G4PVPlacement(0,                       //no rotation
		       G4ThreeVector(xplace,yplace,0),         //at 
		       logicSphere,                //its logical volume
		       materialVec[i],                //its name
		       logicWorld,              //its mother  volume
		       false,                   //no boolean operation
		       0,                       //copy number
		       checkOverlaps);          //overlaps checking
 
     MaterialPosition[i][0] = xplace;
     MaterialPosition[i][1] = yplace;
     MaterialPosition[i][2] = 0.;
  // Set scoring volume to stepping action 
  // (where we will account energy deposit)
  //
     steppingAction->SetVolume(logicSphere);
  }


  // Option to switch on/off checking of volumes overlaps
  //

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
