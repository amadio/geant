
#include "MyDetectorMessenger.hh"

#include "MyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "globals.hh"


MyDetectorMessenger::
MyDetectorMessenger( MyDetectorConstruction* myDet )
  : theDetector( myDet ) {

  theDetectorDir = new G4UIdirectory( "/mydet/" );
  theDetectorDir->SetGuidance( "Detector control." );

  theFieldCommand = new G4UIcmdWithADoubleAndUnit( "/mydet/setField", this );
  theFieldCommand->SetGuidance( "Define uniform magnetic field along Z." );
  theFieldCommand->SetGuidance( " -> in unit of  [Tesla]" );
  theFieldCommand->SetParameterName( "By", false );
  theFieldCommand->SetDefaultValue( 0.0 );
  theFieldCommand->SetUnitCategory( "Magnetic flux density" );
  theFieldCommand->AvailableForStates( G4State_PreInit, G4State_Idle );

  theFieldMapCommand = new G4UIcmdWithABool( "/mydet/useFieldMap", this );
  theFieldMapCommand->SetGuidance( "Define uniform magnetic field along Z." );
  theFieldMapCommand->SetDefaultValue( false );
  theFieldMapCommand->SetParameterName( "UseFieldMap", false );  
  theFieldMapCommand->AvailableForStates( G4State_PreInit, G4State_Idle );
  
  theGDMLCommand = new G4UIcmdWithAString( "/mydet/setGdmlFile", this );
  theGDMLCommand->SetGuidance( "Set the GDML file." );
  theGDMLCommand->SetDefaultValue( "cms2018.gdml" );
  theGDMLCommand->AvailableForStates( G4State_PreInit, G4State_Idle );

}


MyDetectorMessenger::~MyDetectorMessenger() {
  delete theFieldCommand;
  delete theDetectorDir;
}


void MyDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if ( command == theFieldCommand ) {
    theDetector->UseUniformField(true);
    theDetector->SetMagFieldValue( theFieldCommand->GetNewDoubleValue(newValue) );
  }
  if ( command == theFieldMapCommand ) {
    theDetector->UseUniformField( !(theFieldMapCommand->GetNewBoolValue(newValue)) );
  }  
  if ( command == theGDMLCommand ) {
    theDetector->SetGDMLFileName( newValue );
  }

}
