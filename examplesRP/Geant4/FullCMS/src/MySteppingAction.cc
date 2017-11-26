
#include "MySteppingAction.hh"

#include "MyEventAction.hh"
#include "MyTrackInformation.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"


MySteppingAction::MySteppingAction(MyEventAction* evtact) : G4UserSteppingAction(), fEventAction(evtact) {}


MySteppingAction::~MySteppingAction() {}


void MySteppingAction::UserSteppingAction(const G4Step* theStep) {
  G4double edep      = theStep->GetTotalEnergyDeposit();
  G4double stepl     = theStep->GetStepLength();
  G4Track* theTrack  = theStep->GetTrack();
  G4bool   isCharged = (theTrack->GetDefinition()->GetPDGCharge()!=0.);
  G4int primTrackID  = static_cast<MyTrackInformation*>(theTrack->GetUserInformation())->GetPrimaryTrackID();
  fEventAction->AddData(edep, stepl, isCharged, primTrackID-1);
}
