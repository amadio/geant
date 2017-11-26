
#include "MyDetectorMessenger.hh"

#include "MyDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
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

  theGDMLCommand = new G4UIcmdWithAString( "/mydet/setGdmlFile", this );
  theGDMLCommand->SetGuidance( "Set the GDML file." );
  theGDMLCommand->SetDefaultValue( "cms.gdml" );
  theGDMLCommand->AvailableForStates( G4State_PreInit, G4State_Idle );

}


MyDetectorMessenger::~MyDetectorMessenger() {
  delete theFieldCommand;
  delete theDetectorDir;
}


void MyDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if ( command == theFieldCommand ) {
    theDetector->SetMagFieldValue( theFieldCommand->GetNewDoubleValue(newValue) );
  }
  if ( command == theGDMLCommand ) {
    theDetector->SetGDMLFileName( newValue );
  }

}
