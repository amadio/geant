
#ifndef MyEventAction_h
#define MyEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "MyEventDataPerPrimary.hh"

#include <vector>

class G4Event;
class G4Track;

class MyEventAction : public G4UserEventAction {

public:
  MyEventAction();
  virtual ~MyEventAction();

  virtual void BeginOfEventAction(const G4Event *evt);
  virtual void EndOfEventAction(const G4Event *evt);

  void AddData(G4double edep, G4double length, G4bool ischarged, G4int primid);
  void AddSecondaryTrack(const G4Track *track, G4int primid);

private:
  G4int fNumberOfPrimariesPerEvent;
  std::vector<MyEventDataPerPrimary> fEventDataPerPrimary;
};

#endif
