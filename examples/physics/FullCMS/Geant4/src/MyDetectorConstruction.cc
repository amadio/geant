
#include "MyDetectorConstruction.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "MyDetectorMessenger.hh"

G4double MyDetectorConstruction::gFieldValue = 0.0;

MyDetectorConstruction::MyDetectorConstruction()
: fWorld(nullptr), fFieldMgr(nullptr), fUniformMagField(nullptr) , fDetectorMessenger(nullptr) {
  fGDMLFileName = "cms.gdml";
  fFieldValue   = 0.0;
  fDetectorMessenger = new MyDetectorMessenger(this);
}


MyDetectorConstruction::~MyDetectorConstruction() {
  delete fDetectorMessenger;
  if (fUniformMagField) {
    delete fUniformMagField;
  }
}


G4VPhysicalVolume* MyDetectorConstruction::Construct() {
  //  parser.SetOverlapCheck(true);
  fParser.Read(fGDMLFileName,false); // turn off schema checker
  fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fWorld    = (G4VPhysicalVolume *)fParser.GetWorldVolume();
  fWorld->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  if (fWorld==0) {
    G4ExceptionDescription ed;
    ed << "World volume not set properly check your setup selection criteria or GDML input!"
       << G4endl;
    G4Exception( "MyDetectorConstruction::Construct()", "FULL_CMS_0000", FatalException, ed );
  }
  SetMagField();
  return fWorld;
}


void MyDetectorConstruction::SetMagField() {
  if (fUniformMagField ) {
    delete fUniformMagField;
  }
  if (std::abs(fFieldValue)>0.0) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.
    fUniformMagField = new G4UniformMagField(G4ThreeVector(0.0,0.0,fFieldValue));
    fFieldMgr->SetDetectorField(fUniformMagField);
    fFieldMgr->CreateChordFinder(fUniformMagField);
    G4cout << G4endl
           << " *** SETTING MAGNETIC FIELD : fieldValue = " << fFieldValue / tesla
           << " Tesla *** " << G4endl
	         << G4endl;

  } else {
    G4cout << G4endl
           << " *** NO MAGNETIC FIELD SET  *** " << G4endl
	         << G4endl;
  }
}
