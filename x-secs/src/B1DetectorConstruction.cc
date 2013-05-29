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
   char material[98][6] = {"G4_H" ,"G4_He","G4_Li","G4_Be","G4_B" ,"G4_C" ,"G4_N" ,"G4_O" ,"G4_F" ,"G4_Ne",
			  "G4_Na","G4_Mg","G4_Al","G4_Si","G4_P" ,"G4_S" ,"G4_Cl","G4_Ar","G4_K" ,"G4_Ca",
			  "G4_Sc","G4_Ti","G4_V" ,"G4_Cr","G4_Mn","G4_Fe","G4_Co","G4_Ni","G4_Cu","G4_Zn",
			  "G4_Ga","G4_Ge","G4_As","G4_Se","G4_Br","G4_Kr","G4_Rb","G4_Sr","G4_Y" ,"G4_Zr",
			  "G4_Nb","G4_Mo","G4_Tc","G4_Ru","G4_Rh","G4_Pd","G4_Ag","G4_Cd","G4_In","G4_Sn",
			  "G4_Sb","G4_Te","G4_I" ,"G4_Xe","G4_Cs","G4_Ba","G4_La","G4_Ce","G4_Pr","G4_Nd",
			  "G4_Pm","G4_Sm","G4_Eu","G4_Gd","G4_Tb","G4_Dy","G4_Ho","G4_Er","G4_Tm","G4_Yb",
			  "G4_Lu","G4_Hf","G4_Ta","G4_W", "G4_Re","G4_Os","G4_Ir","G4_Pt","G4_Au","G4_Hg",
			  "G4_Tl","G4_Pb","G4_Bi","G4_Po","G4_At","G4_Rn","G4_Fr","G4_Ra","G4_Ac","G4_Th",
			  "G4_Pa","G4_U" ,"G4_Np","G4_Pu","G4_Am","G4_Cm","G4_Bk","G4_Cf"};

   const G4int nmaterials = 98;
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
     G4Material *mat = nist->FindOrBuildMaterial(material[i]);
     G4cout << *mat;
     G4Sphere *solidSphere = new G4Sphere(material[i],             // name
					  0.,kRadius,              // rmin rmax
					  0., 2*pi,                // phi range
					  0., pi);                 // phi range
     G4LogicalVolume* logicSphere =                         
	new G4LogicalVolume(solidSphere,            //its solid
			    mat,                    //its material
			    material[i]);           //its name
  
     G4double xplace = kLattice-kRadius+(i%nballs)*kLattice-0.5*world_sizeXY;
     G4double yplace = kLattice-kRadius+(i/nballs)*kLattice-0.5*world_sizeXY;
     printf("%f %f\n",xplace,yplace);
     new G4PVPlacement(0,                       //no rotation
		       G4ThreeVector(xplace,yplace,0),         //at 
		       logicSphere,                //its logical volume
		       material[i],                //its name
		       logicWorld,              //its mother  volume
		       false,                   //no boolean operation
		       0,                       //copy number
		       checkOverlaps);          //overlaps checking
 
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
