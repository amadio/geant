
#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include "G4GDMLParser.hh"
#include "G4String.hh"

class G4VPhysicalVolume;
class G4FieldManager;
class G4UniformMagField;
class MyDetectorMessenger;

class G4ScalarRZMagFieldFromMap;

class MyDetectorConstruction : public G4VUserDetectorConstruction {

public:
  MyDetectorConstruction();
  ~MyDetectorConstruction();

  G4VPhysicalVolume* Construct() override final;
  void ConstructSDandField() override final;  // needed to create field in each thread   
   // Replaces the old 'SetMagField()'
   
  void SetGDMLFileName(const G4String &gdmlfile) { fGDMLFileName = gdmlfile; }
  void SetMagFieldValue(const G4double fieldValue)
  {
    fFieldValue = fieldValue;
    gFieldValue = fFieldValue;
  }

  static void     UseUniformField(G4bool val= true) {  fUseUniformField= val; }
   
  static G4double GetFieldValue() { return gFieldValue; }
  static G4bool   IsFieldUniform() { return fUseUniformField; }

private:
  // Configuration parameters
  static G4double        gFieldValue;      // primarily used for printout (if uniform...)
  static G4bool          fUseUniformField;

  // 
  G4String                   fGDMLFileName;
  G4double                   fFieldValue;
  G4GDMLParser               fParser;
  G4VPhysicalVolume*         fWorld;
  // G4FieldManager*            fFieldMgr;            // Per-thread !!
  // G4UniformMagField*         fUniformMagField;     // Per-thread !!
  // G4ScalarRZMagFieldFromMap* fSimplifiedCMSfield;  // Per-thread !!
 
  MyDetectorMessenger*       fDetectorMessenger;
};

#endif
