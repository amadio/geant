
#include "MyDetectorConstruction.hh"

#include <atomic>
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4RunManager.hh"
#include "MyDetectorMessenger.hh"
#include "G4ScalarRZMagFieldFromMap.hh"

#include "G4Mag_UsualEqRhs.hh"
#include "G4CashKarpRKF45.hh"
#include "G4ChordFinder.hh"

G4double MyDetectorConstruction::gFieldValue = 0.0;

G4bool   MyDetectorConstruction::fUseUniformField= true;
 
MyDetectorConstruction::MyDetectorConstruction()
<<<<<<< HEAD
: fWorld(nullptr), fFieldMgr(nullptr), fUniformMagField(nullptr),
  fSimplifiedCMSfield(nullptr), fDetectorMessenger(nullptr) {
  fGDMLFileName      = "cms2018.gdml";   
=======
: fWorld(nullptr),
  // fFieldMgr(nullptr), fUniformMagField(nullptr), fSimplifiedCMSfield(nullptr),
  fDetectorMessenger(nullptr) {
  fGDMLFileName = "cms.gdml";
>>>>>>> Enabled choice of field-map or uniform field.  Fixed fieldmap method.
  fFieldValue   = 0.0;
  fDetectorMessenger = new MyDetectorMessenger(this);
}

MyDetectorConstruction::~MyDetectorConstruction()
{
  delete fDetectorMessenger;
<<<<<<< HEAD
}

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
  //  parser.SetOverlapCheck(true);
  fParser.Read(fGDMLFileName, false); // turn off schema checker
  fWorld = (G4VPhysicalVolume *)fParser.GetWorldVolume();
=======
  // delete fUniformMagField;
}

G4VPhysicalVolume* MyDetectorConstruction::Construct() {
  //  parser.SetOverlapCheck(true);
  fParser.Read(fGDMLFileName,false); // turn off schema checker
  // fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fWorld    = (G4VPhysicalVolume *)fParser.GetWorldVolume();
>>>>>>> Enabled choice of field-map or uniform field.  Fixed fieldmap method.
  fWorld->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::Invisible);
  if (fWorld == 0) {
    G4ExceptionDescription ed;
    ed << "World volume not set properly check your setup selection criteria or GDML input!" << G4endl;
    G4Exception("MyDetectorConstruction::Construct()", "FULL_CMS_0000", FatalException, ed);
  }
<<<<<<< HEAD
  // ConstructSDandField();
  return fWorld;
}

// Old method void MyDetectorConstruction::SetMagField() {

void MyDetectorConstruction::ConstructSDandField()
{

  delete fUniformMagField;
  delete fSimplifiedCMSfield;
  fSimplifiedCMSfield= nullptr;
  fUniformMagField = nullptr;
  
#ifndef USE_UNIFORM
  const char *fieldFileName="cmsmagneticfield2015.txt";
  fSimplifiedCMSfield = new G4ScalarRZMagFieldFromMap(fieldFileName);
  fFieldMgr->SetDetectorField(fSimplifiedCMSfield);
  fFieldMgr->CreateChordFinder(fSimplifiedCMSfield);  
#else
  if (std::abs(fFieldValue)>0.0) {
    // Apply a global uniform magnetic field along the Z axis.
    // Notice that only if the magnetic field is not zero, the Geant4
    // transportion in field gets activated.
    auto uniformMagField     = new G4UniformMagField(G4ThreeVector(0.0, 0.0, fFieldValue));
    G4FieldManager *fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
    fieldMgr->SetDetectorField(uniformMagField);
    fieldMgr->CreateChordFinder(uniformMagField);
    G4cout << G4endl << " *** SETTING MAGNETIC FIELD : fieldValue = " << fFieldValue / tesla << " Tesla *** " << G4endl
           << G4endl;

  } else {
    G4cout << G4endl << " *** NO MAGNETIC FIELD SET  *** " << G4endl << G4endl;
=======
  ConstructSDandField();
  return fWorld;
}

static std::atomic<int> callNumConstructField(0);

void MyDetectorConstruction::ConstructSDandField() {
  G4MagneticField* magField= nullptr;
  const char *fieldFileName="cmsmagfield2015.txt";  
  auto fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
   
  if( fUseUniformField ) {
    if (std::abs(fFieldValue)>0.0) {
      // Apply a global uniform magnetic field along the Z axis.
      // Notice that only if the magnetic field is not zero, the Geant4
      // transportion in field gets activated.
      auto uniformMagField = new G4UniformMagField(G4ThreeVector(0.0,0.0,fFieldValue));
      magField = uniformMagField;
      
      G4cout << G4endl
             << " *** SETTING MAGNETIC FIELD : fieldValue = " << fFieldValue / tesla
             << " Tesla *** " << G4endl
             << G4endl;
    } else {
       G4cout << G4endl
              << " *** NO MAGNETIC FIELD SET  *** " << G4endl
              << G4endl;
    }
  } else {
    auto simplifiedCMSfield = new G4ScalarRZMagFieldFromMap(fieldFileName);
    magField = simplifiedCMSfield;

    double  position[4] = { 0., 0., 0., 0. }, fieldArr[3];
    simplifiedCMSfield->GetFieldValue( position, fieldArr );

    // for( int i=0; i<3; i++) { fieldArr[i] *= 1.0 / tesla; }
    G4ThreeVector fieldVec( fieldArr[0], fieldArr[1], fieldArr[2] );
    SetMagFieldValue( fieldArr[2] ); // Know that it is solinoidal - B // Bz at origin
    fieldVec *= 1.0 / tesla;
    G4cout << G4endl
           << " *** SET MAGNETIC FIELD : Scalar RZ Mag Field From Map.  Value at origin = "
        // << fieldArr[0] << " , "  << fieldArr[1] << " , "  << fieldArr[2]
           << fieldVec.x() << " , "  << fieldVec.y() << " , "  << fieldVec.z()       
           << " Tesla *** " << G4endl;
>>>>>>> Enabled choice of field-map or uniform field.  Fixed fieldmap method.
  }
  // G4cout << " Dbg:  fieldMgr = " << fieldMgr << G4endl;

  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);
  /*****
  auto equation = new G4Mag_UsualEqRhs(magField);
  auto Stepper = new G4CashKarpRKF45( equation );
  G4cout<<"  - Using G4CashKarpRKF45 stepper. " << G4endl;

  double minStep= 1.0e-2 * CLHEP::millimeter;  // Geant4 default - to be configured
  
  auto chordFinder = new G4ChordFinder( magField, minStep, Stepper );
  fieldMgr->SetChordFinder( chordFinder );
  G4cout << " - Created and set ChordFinder. " << G4endl;
   *****/
  G4cout << " ConstructSDandField: Dbg>  Call Number " << callNumConstructField++ << G4endl;
}
