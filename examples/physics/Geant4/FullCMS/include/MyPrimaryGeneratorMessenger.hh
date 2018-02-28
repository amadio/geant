
#ifndef MyPrimaryGeneratorMessenger_h
#define MyPrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class MyPrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;

class MyPrimaryGeneratorMessenger : public G4UImessenger {
public:
  MyPrimaryGeneratorMessenger(MyPrimaryGeneratorAction *gun);
  ~MyPrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand *, G4String);

private:
  MyPrimaryGeneratorAction *fTheGun;

  G4UIdirectory *fGunDirectory;
  G4UIcmdWithAnInteger *fNumPrimaryPerEvtCmd;
  G4UIcmdWithAString *fPrimaryTypeCmd;
  G4UIcmdWithADoubleAndUnit *fPrimaryEnergyCmd;
  G4UIcmdWith3Vector *fPrimaryDirCmd;
};

#endif
