
#include "MyTrackingAction.hh"

#include "G4Track.hh"
#include "G4TrackVector.hh"
#include "G4TrackingManager.hh"
#include "MyEventAction.hh"
#include "MyTrackInformation.hh"

MyTrackingAction::MyTrackingAction(MyEventAction *evtact) : G4UserTrackingAction(), fEventAction(evtact)
{
}

MyTrackingAction::~MyTrackingAction()
{
}

void MyTrackingAction::PreUserTrackingAction(const G4Track *theTrack)
{
  if (theTrack->GetParentID() == 0) { // primary track so set its auxiliary track member: primary trackID
    theTrack->SetUserInformation(new MyTrackInformation(theTrack->GetTrackID()));
  } else {
    G4int primTrackID = static_cast<MyTrackInformation *>(theTrack->GetUserInformation())->GetPrimaryTrackID();
    fEventAction->AddSecondaryTrack(theTrack, primTrackID - 1);
  }
}

void MyTrackingAction::PostUserTrackingAction(const G4Track *theTrack)
{
  G4TrackVector *secondaries = fpTrackingManager->GimmeSecondaries();
  if (secondaries) {
    MyTrackInformation *info = static_cast<MyTrackInformation *>(theTrack->GetUserInformation());
    G4int primaryTrackID     = info->GetPrimaryTrackID(); // get the primary trackID
    G4int currentTrackID     = theTrack->GetTrackID();    // track ID of theTrack
    size_t nSecondaries      = secondaries->size();
    for (size_t isec = 0; isec < nSecondaries; ++isec) {
      // if the secondary was created by theTrack set its auxiliary info to store the primaryTrackID
      G4Track *secTrack = (*secondaries)[isec];
      if (secTrack->GetParentID() == currentTrackID) {
        secTrack->SetUserInformation(new MyTrackInformation(primaryTrackID));
      }
    }
  }
}
