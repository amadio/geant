
#include "MyPrimaryGeneratorMessenger.hh"

#include "MyPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"


MyPrimaryGeneratorMessenger::MyPrimaryGeneratorMessenger(MyPrimaryGeneratorAction *gun) : G4UImessenger(),
 fTheGun(gun), fGunDirectory(nullptr), fPrimaryTypeCmd(nullptr), fPrimaryEnergyCmd(nullptr), fPrimaryDirCmd(nullptr) {
 fGunDirectory     = new G4UIdirectory("/mygun/");
 fGunDirectory->SetGuidance("gun control");

 fNumPrimaryPerEvtCmd = new G4UIcmdWithAnInteger("/mygun/primaryPerEvt",this);
 fNumPrimaryPerEvtCmd->SetGuidance("set number of primary particles per event");
 fNumPrimaryPerEvtCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 fPrimaryTypeCmd      = new G4UIcmdWithAString("/mygun/particle",this);
 fPrimaryTypeCmd->SetGuidance("set primary particle type");
 fPrimaryTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 fPrimaryEnergyCmd    = new G4UIcmdWithADoubleAndUnit("/mygun/energy",this);
 fPrimaryEnergyCmd->SetGuidance("set primary particle energy");
 fPrimaryEnergyCmd->SetParameterName("pEnergy",false);
 fPrimaryEnergyCmd->SetRange("pEnergy>0.");
 fPrimaryEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 fPrimaryDirCmd       = new G4UIcmdWith3Vector("/mygun/direction",this);
 fPrimaryDirCmd->SetGuidance("set primary particle direction");
 fPrimaryDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}


MyPrimaryGeneratorMessenger::~MyPrimaryGeneratorMessenger() {
  delete fGunDirectory;
  delete fNumPrimaryPerEvtCmd;
  delete fPrimaryTypeCmd;
  delete fPrimaryEnergyCmd;
  delete fPrimaryDirCmd;
}

void MyPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command==fNumPrimaryPerEvtCmd) {
    fTheGun->SetNumPrimaryPerEvt(fNumPrimaryPerEvtCmd->GetNewIntValue(newValue));
  }
  if (command==fPrimaryTypeCmd) {
    fTheGun->SetPrimaryName(newValue);
  }
  if (command==fPrimaryEnergyCmd) {
    fTheGun->SetPrimaryEnergy(fPrimaryEnergyCmd->GetNewDoubleValue(newValue));
  }
  if (command==fPrimaryDirCmd) {
    fTheGun->SetPrimaryDirection(fPrimaryDirCmd->GetNew3VectorValue(newValue));
  }
}
