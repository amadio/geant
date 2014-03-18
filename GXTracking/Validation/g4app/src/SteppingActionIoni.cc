//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingActionIoni.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4ThreeVector.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4Material.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"

#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4GammaConversion.hh"

#include "G4dataBuffer.hh"
#include "GPHistoManager.hh"
#include <vector>
using std::vector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingActionIoni::SteppingActionIoni()                                         
{
  detector = (DetectorConstruction*)
             G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  eventaction = (EventAction*)
                G4RunManager::GetRunManager()->GetUserEventAction();               

  first_eBrem = true;
  first_eIoni = true;
  first_eMsc = true;

  first_gCompt = true;
  first_gPhot = true;
  first_gConv = true;

 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingActionIoni::~SteppingActionIoni()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingActionIoni::UserSteppingAction(const G4Step* aStep)
{
  // get volume of the current step
  //G4VPhysicalVolume* volume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // collect energy and track length step by step
  //G4double edep = aStep->GetTotalEnergyDeposit();
  
  G4double stepl = 0.;
  G4Track* track = aStep->GetTrack();
  if(track->GetDefinition()->GetPDGCharge() != 0.) stepl = aStep->GetStepLength();
  G4double mass = track->GetDefinition()->GetPDGMass();
  
  if(/*storeTable*/ 1) {
    const G4Material* mat = aStep->GetTrack()->GetMaterial();
    //  write physics tables only for the crystal
    if( mat->GetName() == "PbWO4") WritePhysicsTables(aStep);
  }

  // kill every track, so each  track takes a single step
  track->SetTrackStatus(fStopAndKill);

  // G4cout <<"SteppingActionIoni: trkID="<< track->GetTrackID()
  // 	 <<" energy="<< track->GetTotalEnergy()
  // 	 <<" pos="<< track->GetPosition()
  // 	 <<" momDirection="<< track->GetMomentumDirection() << G4endl;

  //=============================================
#ifdef GPUPLOTS
  GPHistoManager& hmgr = GPHistoManager::getInstance();
#endif

  //== make some plots: postStep electron
  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4double eleEne = preStepPoint->GetTotalEnergy();
  G4double elemom = sqrt(eleEne*eleEne - mass*mass);
  G4ThreeVector momVec = elemom*preStepPoint->GetMomentumDirection();
#ifdef GPUDEBUG
  G4cout<<"\nBefore Ioni: Ene="<< eleEne/MeV <<" Ekin="<< (eleEne-mass)/MeV <<' '<< preStepPoint->GetKineticEnergy()/MeV
	<<" momVec=("<< momVec.x()/MeV <<"; "<< momVec.y()/MeV <<"; "<< momVec.z()/MeV <<")" << G4endl;
#endif
  G4double elepx = momVec.x();
  G4double elepy = momVec.y();
  G4double elepz = momVec.z();
  G4double inputMomentum = sqrt(elepx*elepx + elepy*elepy + elepz*elepz);

  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();
  eleEne = postStepPoint->GetTotalEnergy();
  elemom = sqrt(eleEne*eleEne - mass*mass);
  momVec = elemom*postStepPoint->GetMomentumDirection();

  //== make some plots: electron change due to physics process
  G4double pcEtot = aStep->GetPostStepPoint()->GetTotalEnergy();
  const G4ThreeVector& pcdir = aStep->GetPostStepPoint()->GetMomentumDirection();
  G4double pcdirx = pcdir.x();
  G4double pcdiry = pcdir.y();
  G4double pcdirz = pcdir.z();
  G4double pcAngle = acos((elepx*pcdirx+elepy*pcdiry+elepz*pcdirz)/inputMomentum);

#ifdef GPUPLOTS
  hmgr.getHisto("hioniEnergy").Fill(pcEtot-mass);
  hmgr.getHisto("hioniAngle").Fill(pcAngle);

  G4dataBuffer& buffers = G4dataBuffer::getInstance();
#ifdef GPUDEBUG
  G4cout<<"SteppingActionIoni: preStepLambda buffer: "<< buffers.getBuffer("hioniPreStepLambda").size() << G4endl;
#endif

  hmgr.getHisto("hioniStepLengthPost").Fill(stepl/mm);
  hmgr.getHisto("hioniEnergyLoss").Fill( (preStepPoint->GetTotalEnergy()-postStepPoint->GetTotalEnergy())/MeV );

  hmgr.fillFromVector("hioniPreStepLambda", buffers.getBuffer("hioniPreStepLambda") );
  hmgr.fillFromVector("hioniPreStepScaledEnergy", buffers.getBuffer("hioniPreStepScaledEnergy") );
  hmgr.fillFromVector("hioniNbOfIntLengthLeft", buffers.getBuffer("hioniNbOfIntLengthLeft") );

  hmgr.fillFromVector("hioniDedxForScaledEnergyTimesLength", buffers.getBuffer("hioniDedxForScaledEnergyTimesLength"));
  hmgr.fillFromVector("hioniElossFromKinEnergyMinusScaledEnergyForLoss", buffers.getBuffer("hioniElossFromKinEnergyMinusScaledEnergyForLoss"));
  hmgr.fillFromVector("hioniElossFromSampleFluctuations", buffers.getBuffer("hioniElossFromSampleFluctuations"));
  hmgr.fillFromVector("hioniEloss", buffers.getBuffer("hioniEloss"));
#endif // GPUPLOTS

#ifdef GPUDEBUG
  G4cout<<"After ioni: Ene="<< eleEne/MeV <<" Ekin="<< (eleEne-mass)/MeV <<' '<< postStepPoint->GetKineticEnergy()/MeV
	<<" momVec=("<< momVec.x()/MeV <<"; "<< momVec.y()/MeV <<"; "<< momVec.z()/MeV <<")"
	<<" pcEtot="<< pcEtot << " pcAngle="<< pcAngle << G4endl;
#endif

  // G4double pcMom = sqrt(pcEtot*pcEtot - mass*mass);
  // printf(" outElectron: P[MeV] = (%e;%e;%e)\n",pcMom*pcdirx,pcMom*pcdiry,pcMom*pcdirz);

  // secondaries
  const vector<const G4Track*>* secTracks = aStep->GetSecondaryInCurrentStep();
  vector<const G4Track*>::const_iterator itrack = secTracks->begin();

  for( ; itrack!=secTracks->end(); ++itrack) {
    //offset is a global counter for the last array position of secondaries 
    const G4Track* sectrk = *itrack;
    // G4cout<<"   Secondary: trkID="<< sectrk->GetTrackID()
    //   <<" energy="<< sectrk->GetTotalEnergy()
    //   <<" pos="<< sectrk->GetPosition()
    //   <<" momDirection="<< sectrk->GetMomentumDirection() << G4endl;

    const G4ThreeVector& secp = sectrk->GetMomentum();
    G4double secpx = secp.x();
    G4double secpy = secp.y();
    G4double secpz = secp.z();
    G4double secQ = sectrk->GetParticleDefinition()->GetPDGCharge();
    G4double secmass = mass*secQ*secQ;
    G4double secmom = sqrt(secpx*secpx+secpy*secpy+secpz*secpz);
    G4double secKinE = sqrt(secmom*secmom+secmass*secmass)-secmass;
    G4double angle = acos((elepx*secpx+elepy*secpy+elepz*secpz)/(inputMomentum*secmom));

#ifdef GPUDEBUG
    printf("Gamma: E=%e,  P=(%e;%e;%e)\n",secmom/MeV,secpx/MeV,secpy/MeV,secpz/MeV);
    printf("Mom.Cons.uncertainty: P=(%e;%e;%e)\n",(elepx-elemom*pcdirx-secpx)/elepx*100,
     	   (elepy-elemom*pcdiry-secpy)/elepy*100,(elepz-elemom*pcdirz-secpz)/elepz*100);
    printf("Angles: pc=%e, gamma=%e\n\n",pcAngle,angle);
#endif

#ifdef GPUPLOTS
    hmgr.getHisto("hioniSecAngle").Fill(angle);
    hmgr.getHisto("hioniSecEnergy").Fill(secKinE/MeV);
#endif
  }
  //=============================================

  // if (volume == detector->GetAbsorber()) eventaction->AddAbs(edep,stepl);
  // if (volume == detector->GetGap())      eventaction->AddGap(edep,stepl);
  
  //example of saving random number seed of this event, under condition
  //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingActionIoni::WritePhysicsTables( const G4Step * theStep ) {

  const G4ParticleDefinition* part = theStep->GetTrack()->GetDefinition();

  G4ProcessManager* pm = 
    theStep->GetTrack()->GetDefinition()->GetProcessManager();
  G4ProcessVector* pv = pm->GetProcessList();

  int np = pv->length();

  for(int i = 0; i < np ; ++i) {

    G4VProcess* proc = (*pv)[i];
    const G4String& name = proc->GetProcessName();

    //electrons
    if(name == "eBrem" && first_eBrem ) {
      G4eBremsstrahlung* eproc = (G4eBremsstrahlung*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_eBrem = false;
    }
    if(name == "eIoni" && first_eIoni) {
      G4eIonisation* eproc = (G4eIonisation*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_eIoni = false;
    }
    if(name == "msc" && first_eMsc) {
      G4eMultipleScattering* eproc = (G4eMultipleScattering*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_eMsc = false;
    }

    //photons
    if(name == "compt" && first_gCompt) {
      G4ComptonScattering* eproc = (G4ComptonScattering*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_gCompt = false;
    }
    if(name == "phot" && first_gPhot) {
      G4PhotoElectricEffect* eproc = (G4PhotoElectricEffect*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_gPhot = false;
    }
    if(name == "conv" && first_gConv) {
      G4GammaConversion* eproc = (G4GammaConversion*) proc;
      G4bool yes = eproc->StorePhysicsTable(part,"./table",true);
      if(yes) first_gConv = false;
    }
  }
}
