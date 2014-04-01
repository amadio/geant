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

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class cmsExpMagneticField;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction {

public:
  
  DetectorConstruction(G4int nphi, G4int nz, G4double density);
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

public:
     
  void SetupMaterial (G4String);     
  void SetDimensions(G4double);

  void SetMagField(G4double);
     
  void UpdateGeometry();
  //void PrintDetectorInfo(); 

  // G4double GetWorldSizeXY()    {return _worldSizeX;}; 
  // G4double GetWorldSizeZ()     {return _worldSizeYZ;};
  // G4double GetCalorThickness() {return CalorThickness;};
  // G4double GetCalorSizeZ()     {return CalorSizeZ;};

  //G4int GetNLayers()           {return _nLayers;};
  const G4VPhysicalVolume* GetPhysWorld() {return _physWorld;};           

private:

  G4int   _nPhi;
  G4int   _nZ;
  G4double _density;

  G4Material*        _defaultMaterial;
  G4Material*        _crystalMaterial;

  G4Box*             _solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   _logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* _physWorld;    //pointer to the physical World

  G4Tubs**            _solidVols;
  G4LogicalVolume**   _logicVols;
  G4VPhysicalVolume** _physVols;

  //  G4UniformMagField* _magField;      //pointer to the magnetic field
  cmsExpMagneticField* _magField;      //pointer to the magnetic field
  DetectorMessenger* _detectorMessenger;  //pointer to the Messenger

private:

  void DefineMaterials();
  G4VPhysicalVolume* ConstructCalorimeter();
};

#endif

