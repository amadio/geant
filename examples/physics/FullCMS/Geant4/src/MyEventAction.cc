
#include "MyEventAction.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4ThreeVector.hh"
#include "G4PrimaryVertex.hh"

#include "MyRunAction.hh"
#include "MySteppingAction.hh"
#include "MyRun.hh"
#include "MyPrimaryGeneratorAction.hh"


MyEventAction::MyEventAction() : G4UserEventAction() {
  fNumberOfPrimariesPerEvent = 10; // will be changed automatically if needed
  fEventDataPerPrimary.resize(fNumberOfPrimariesPerEvent);
  for (G4int id=0; id<fNumberOfPrimariesPerEvent; ++id) {
    fEventDataPerPrimary[id].Clear();
  }
}


MyEventAction::~MyEventAction() {
  fEventDataPerPrimary.clear();
}


void MyEventAction::BeginOfEventAction(const G4Event* evt) {
  G4int curNumPrims = evt->GetNumberOfPrimaryVertex();
  if (curNumPrims>fNumberOfPrimariesPerEvent) {
    fNumberOfPrimariesPerEvent = curNumPrims;
    fEventDataPerPrimary.resize(fNumberOfPrimariesPerEvent);
  }
  for (G4int id=0; id<curNumPrims; ++id) {
    fEventDataPerPrimary[id].Clear();
  }
}


void MyEventAction::EndOfEventAction(const G4Event* evt) {
  //get the Run and add the data collected during the simulation of the event that has been completed
  MyRun* run = static_cast<MyRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4int  nPrimPart = evt->GetNumberOfPrimaryVertex();
  G4cout << " \n================================================================= \n"
         << " ===  EndOfEventAction  --- event = " << evt->GetEventID() << " with "<< nPrimPart << " primary:"<<G4endl;
  for (G4int ip=0; ip<nPrimPart; ++ip) {
    G4PrimaryParticle*   pPart      = evt->GetPrimaryVertex(ip)->GetPrimary(0);
    G4String             pPartName  = pPart->GetParticleDefinition()->GetParticleName();
    G4double             pPartEkin  = pPart->GetKineticEnergy();
    const G4ThreeVector& pPartDir   = pPart->GetMomentumDirection();
    G4int primaryTypeIndx = MyPrimaryGeneratorAction::GetPrimaryTypeIndex(pPartName);
    run->FillPerEventNumPrimaries(primaryTypeIndx);
    run->FillPerEvent(fEventDataPerPrimary[ip],primaryTypeIndx);
    G4cout<< "  Primary Particle:  " << ip  << " (type inedx = " << primaryTypeIndx << ")\n"
          << "    Name      =  "     << pPartName                                   << " \n"
          << "    Energy    =  "     << pPartEkin/GeV                          << " [GeV]\n"
          << "    Direction =  "     << pPartDir                                    << " \n"
          << fEventDataPerPrimary[ip]<< G4endl;
  }
}


void MyEventAction::AddData(G4double edep, G4double length, G4bool ischarged, G4int primid) {
  fEventDataPerPrimary[primid].fEdep += edep;
  if (ischarged) {
    fEventDataPerPrimary[primid].fTrackLCh    += length;
    fEventDataPerPrimary[primid].fChargedStep ++;
  } else {
    fEventDataPerPrimary[primid].fTrackLNe    += length;
    fEventDataPerPrimary[primid].fNeutralStep ++;
  }
}


void MyEventAction::AddSecondaryTrack(const G4Track* track, G4int primid) {
  const G4ParticleDefinition* pdf = track->GetDefinition();
  if (pdf==G4Gamma::Gamma()) {
    fEventDataPerPrimary[primid].fNGamma += 1.;
  } else if (pdf==G4Electron::Electron()) {
    fEventDataPerPrimary[primid].fNElec += 1.;
  } else if (pdf==G4Positron::Positron()) {
    fEventDataPerPrimary[primid].fNPosit += 1.;
  }
}


void MyEventAction::AddGlobalTime(G4double gtime, bool ischarged, G4int primid) {
  if (ischarged) {
    fEventDataPerPrimary[primid].fTrackTCh += gtime;
  } else {
    fEventDataPerPrimary[primid].fTrackTNe += gtime;
  }
}
