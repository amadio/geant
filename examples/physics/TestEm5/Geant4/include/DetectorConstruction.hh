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
/// \file electromagnetic/TestEm5/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 98752 2016-08-09 13:44:40Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class DetectorMessenger;

class G4UniformMagField;
class G4FieldManager;
class PrimaryGeneratorAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:
 
     virtual G4VPhysicalVolume* Construct();

     void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);
     void SetAbsorberSizeYZ   (G4double);

     void SetAbsorberXpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeX   (G4double);
     void SetWorldSizeYZ  (G4double);

     void SetMagField(const G4ThreeVector& fv) { fMagFieldVector = fv; }

     void SetPrimaryGenerator(PrimaryGeneratorAction* pg) { fPrimaryGenerator = pg; }

  public:

     void PrintCalorParameters();

     G4Material* GetAbsorberMaterial()  {return fAbsorberMaterial;};
     G4double    GetAbsorberThickness() {return fAbsorberThickness;};
     G4double    GetAbsorberSizeYZ()    {return fAbsorberSizeYZ;};

     G4double    GetAbsorberXpos()      {return fXposAbs;};
     G4double    GetxstartAbs()         {return fXstartAbs;};
     G4double    GetxendAbs()           {return fXendAbs;};

     G4Material* GetWorldMaterial()     {return fWorldMaterial;};
     G4double    GetWorldSizeX()        {return fWorldSizeX;};

     const G4VPhysicalVolume* GetAbsorber() {return fPhysiAbsorber;};

  private:

     void ChangeGeometry();

     G4Material*        fAbsorberMaterial;
     G4double           fAbsorberThickness;
     G4double           fAbsorberSizeYZ;

     G4double           fXposAbs;
     G4double           fXstartAbs, fXendAbs;

     G4Material*        fWorldMaterial;
     G4double           fWorldSizeX;
     G4double           fWorldSizeYZ;

     G4bool             fDefaultWorld;

     G4Box*             fSolidWorld;
     G4LogicalVolume*   fLogicWorld;
     G4VPhysicalVolume* fPhysiWorld;

     G4Box*             fSolidAbsorber;
     G4LogicalVolume*   fLogicAbsorber;
     G4VPhysicalVolume* fPhysiAbsorber;
     
     // field related members
     G4ThreeVector      fMagFieldVector;
     G4FieldManager*    fFieldMgr;
     G4UniformMagField* fUniformMagField;

     PrimaryGeneratorAction* fPrimaryGenerator;
     
     DetectorMessenger* fDetectorMessenger;

  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     void SetConstantField();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

