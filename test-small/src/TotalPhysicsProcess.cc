//
// Authors: M. Nowak, S.Y. Jun & J. Apostolakis, April 2014
//
//
// Started from TabulatedPhysicsProcess

#include "TotalPhysicsProcess.hh"
#include "TabulatedDataManager.hh"
#include "TabulatedPhysicsNameSpace.hh"

#if 0
#include "TTabPhysMgr.hh"
TotalPhysicsProcess::TTabPhysMgr* theTabPhysManager= 0;
#endif

#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4GPILSelection.hh"
#include "G4SystemOfUnits.hh"

#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

/*
#include "TSystem.h"
#include <TPartIndex.h>
*/

#include "MaterialConverter.hh"

#include "TParticlePDG.h"
#include "TDatabasePDG.h"
#include "TPartIndex.h"

#include "TEXsec.h"
#include "TEFstate.h"

// TotalPhysicsProcess::TabulatedDataManager* theDataManager= 0;

// 3 keV cut in energy by default
G4double TotalPhysicsProcess::fgEnergyLimit = 3.0e-6;  

TotalPhysicsProcess::TotalPhysicsProcess(G4String processName) :
  G4VRestContinuousDiscreteProcess(processName,fGeneral),
  fMaterialIndex(-1)
  // fParticleChange(0)
{
  // fSecDefinition = G4Geantino::Geantino();
  // theDataManager = 0;
  theDataManager = TabulatedDataManager::Instance();
  if( ! theDataManager ){
    G4cerr << "ERROR TotalPhysicsProcess constructor> Cannot obtain instance of Tabulated DAta Manager. Fatal Error." << G4endl;
    static int count =0;
    if ( ++ count > 5 ) { exit(1); }
  }
  fSecParticles.reserve(5);
  fParticleChange= new G4ParticleChange();
}

TotalPhysicsProcess::~TotalPhysicsProcess() {
  delete fParticleChange;
}

int TotalPhysicsProcess::SetupForMaterial(const G4Track& track)
{
  //get the material index for the current volume

  G4Material* materialG4= track.GetMaterial();

  static MaterialConverter *sMatConverter= MaterialConverter::Instance(); 
  fMaterialIndex= sMatConverter->GetRootMaterial( materialG4->GetIndex() ) ;
 
  return fMaterialIndex; 
}

// #include "TGeoRCExtension.h"
#include "TMXsec.h"
#include "TEXsec.h"
#include  "TList.h"

//______________________________________________________________________________ 
G4double TotalPhysicsProcess::GetContinuousStepLimit(const G4Track& track, 
                                                     G4double previousStepSize,
                                                     G4double currentMinimumStep,
                                                     G4double &currentSafety){
    G4int rootMatId = SetupForMaterial(track);
    theDataManager->ApplyMsc(rootMatId, track);
    return DBL_MAX;
}

//_____________________________________________________________________________
G4VParticleChange* TotalPhysicsProcess::AlongStepDoIt(const G4Track& track, 
                                                      const G4Step& step){
   // init G4ParticleChange from current track
  fParticleChange->Initialize(track);
  // get rootgeom corresponding material index 
  G4int rootMatId= SetupForMaterial(track);
  // compute energy loss from dedx and update the particle change; set status to
  // fAlive, fStopButAlive or fStopAndKill if final Ekin> fgEnergyLimit, <fgEnergyLimit
  // and it doesn't have or has NuclearCaptureAtRest   
  theDataManager->EnergyLoss(rootMatId, track, step, fParticleChange, fgEnergyLimit);   

  return fParticleChange;
}

//______________________________________________________________________________
G4VParticleChange* TotalPhysicsProcess::AtRestDoIt(const G4Track &atrack, 
                                                    const G4Step &astep){
   // DOES NOTHING AT THE MOMENT JUST KILLS THE TRACK!
   // init G4ParticleChange from current track
   fParticleChange->Initialize(atrack); 

   // get rootgeom corresponding material index 
   G4int rootMatId = SetupForMaterial(atrack);

   // Sample final state for at rest process (we have only uclear capture at 
   // rest) and fill the particle change 
   theDataManager->SampleFinalStateAtRest(rootMatId, atrack, fParticleChange, 
                                          fgEnergyLimit);
   return fParticleChange;
}

//______________________________________________________________________________
G4double TotalPhysicsProcess::GetMeanFreePath(const G4Track& track,
                                              G4double , // previousStepSize,
                                              G4ForceCondition* condition){
  // get rootgeom corresponding material index 
  G4int rootMatId= SetupForMaterial(track);

  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double preStepLambda =
     theDataManager->GetInteractionLength(rootMatId,track); 
  return preStepLambda;
}

// Invoke post-step-do-it action
//______________________________________________________________________________
G4VParticleChange* TotalPhysicsProcess::PostStepDoIt(const G4Track& track, 
                                                     const G4Step& step) {

  // init G4ParticleChange from current track
  fParticleChange->Initialize(track);

  // Check if it has survived the along-step-do-it part:
  // If yes: early return
  // If no : go on
  if( track.GetTrackStatus() != fAlive )
    return fParticleChange;


  // Sampling element for interaction and type of interaction on that
  Int_t reactionId   = -1;
  Int_t elementIndex = -1;

  G4int rootMatId= SetupForMaterial(track);
  elementIndex = theDataManager->SampleInteraction(rootMatId, track, 
                                                   reactionId);

  // Go for the corresponding final states: THIS IS THE POINT WHERE TABULATED
  // PHYSCICS DATA CAN BE REPLACED BY VECTORIZED FULL DISCRETE PHYSCICS LATER!!!
  if( reactionId < 0 ) { // if there is no reaction for this particle
    fParticleChange->SetNumberOfSecondaries(0);
    // update primary information; not realy necessary since it was initialize
    fParticleChange->ProposeLocalEnergyDeposit(0.);
    fParticleChange->ProposeNonIonizingEnergyDeposit(0.);
    fParticleChange->ProposeTrackStatus(fAlive);
    return fParticleChange; 
  }    

  // Sample final state for the sampled interaction and fill the particle change 
  theDataManager->SampleFinalState(elementIndex, reactionId, track, 
                                   fParticleChange, fgEnergyLimit );
  
  return fParticleChange;
}




//This is not for MISI
//______________________________________________________________________________
void TotalPhysicsProcess::Print(const G4Step& step) 
{
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  std::cout << "TotalPhysicsProcess ProcName, E, dE, x, StepLength Nsec " 
	    << postStepPoint->GetProcessDefinedStep()->GetProcessName() << " " 
	    << postStepPoint->GetKineticEnergy()/MeV  << " " 
	    << step.GetTotalEnergyDeposit()/MeV  << " " 
	    << postStepPoint->GetPosition() << " "  
	    << step.GetStepLength() << " " 
	    << fParticleChange->GetNumberOfSecondaries() << std::endl;
}

#include "G4ParticleTable.hh"

const G4ParticleDefinition* TotalPhysicsProcess::ParticleDefinition(G4int ipdg) 
{
  //only for gamma/electron/photon for the standard EM processes
  const G4ParticleDefinition* pdef;
  static G4ParticleTable* theParticleTable = G4ParticleTable::GetParticleTable();
  
  switch (ipdg) {
  case 22:
    pdef = G4Gamma::Gamma(); 
    break;
  case 11:
    pdef = G4Electron::Electron(); 
    break;
  case -11:
    pdef = G4Positron::Positron(); 
    break;
  default:
      pdef = theParticleTable->FindParticle(ipdg);
    break;
  }

  return pdef;
}
