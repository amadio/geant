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



TotalPhysicsProcess::TotalPhysicsProcess(G4String processName) :
  G4VDiscreteProcess(processName,fGeneral),
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

// G4double TotalPhysicsProcess::MeanFreePath(const G4Track& track)

G4double TotalPhysicsProcess::GetMeanFreePath(const G4Track& track,
                                              G4double , // previousStepSize,
                                              G4ForceCondition* condition) 
{
  // G4double energy = track.GetKineticEnergy()/GeV;

  int rootMatId= SetupForMaterial(track);

  //Find the number of particles with reactions
  //
  int pdgEncoding= track.GetParticleDefinition()->GetPDGEncoding();

  // Find the equivalent Root particle
  static TPartIndex* partIdx= TPartIndex::I();
  fParticleId = partIdx->PartIndex( pdgEncoding ); // GeantVparticle index
  
  // std::cout << " Particle Indices:  pdg= " << pdgEncoding << " rootId = "  << partId << std::endl;
  // std::cout << " Number of Particles: " << partIdx->NPart() << std::endl;
  
  // condition is set to "Not Forced"
  *condition = NotForced;
  
  G4double preStepLambda =
     theDataManager->GetInteractionLength(rootMatId,track); // partId,energy);
  return preStepLambda;
}

G4VParticleChange* 
TotalPhysicsProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{

  // Sampling element for interaction and type of interaction on that
  Int_t reactionId   = -1;
  Int_t elementIndex = -1;

  elementIndex = theDataManager->SampleInteraction(fMaterialIndex, track, 
                                                   reactionId);

  // Go for the corresponding final states: THIS IS THE POINT WHERE TABULATED
  // PHYSCICS DATA CAN BE REPLACED BY VECTORIZED FULL DISCRETE PHYSCICS LATER!!!
  if( reactionId < 0 ) { // if there is no reaction for this particle
    fParticleChange->SetNumberOfSecondaries(0);
    // update primary information
    fParticleChange->ProposeLocalEnergyDeposit(0.);
    fParticleChange->ProposeNonIonizingEnergyDeposit(0.);
    fParticleChange->ProposeTrackStatus(fAlive);
 
    // TO DO:: need to be sure that everything remain the same
   

    fParticleChange->ProposeEnergy(track.GetKineticEnergy());
    fParticleChange->ProposeMomentumDirection(track.GetMomentumDirection());
    fParticleChange->ProposeProperTime(track.GetProperTime());
    fParticleChange->ProposePosition(track.GetPosition());
    fParticleChange->ProposeGlobalTime(track.GetGlobalTime());

    fParticleChange->ProposeMass(track.GetDynamicParticle()->GetMass());
    fParticleChange->ProposeCharge(track.GetDynamicParticle()->GetCharge());

    return fParticleChange; 
  }    

  theDataManager->SampleFinalState(elementIndex, reactionId, track, 
                                   fParticleChange);
  
  return fParticleChange;

/* 
  I PUT THE WHOLE EXISTING IMPLEMENTATION INTO COMMENT BECAUSE THERE ARE MANY 
  PROBLEMS HERE RANGING FROM SIMPLE ZERO POINTER EXCEPTION TO SEMANTICALLY 
  INCORRECT, INCOMPLETE IMPLEMENTATION. SEE THE CORRECT IMPLEMENTATION ABOVE AND 
  THE CORRECT IMPLEMENTATION OF THE NECESSARY METHODS IN TabulatedDataManager!!!

  // Int_t TTabPhysMgr::SampleInt(Int_t imat, Int_t ntracks, GeantTrack_v &tracks, Int_t tid)
  
  // static 
  TGeoManager* tGeom= MaterialConverter::GetTGeomManager();
  fSecParticles.clear();

  // OLD, to be replaced - see below
  // G5proc   fReaction= kTotal;
  theDataManager->SampleSecondaries(&fSecParticles,fMaterialIndex,
				    &track,kTotal);
  
  // TODO:  Replace the following method with one which samples ALL reactions - not just the one fReaction
  
  TGeoMaterial *mat;
  mat = (TGeoMaterial *) (tGeom->GetListOfMaterials())->At(fMaterialIndex);
  
  // Find the relevant TMXsec object for this particle, material (?)
  // TMXsec = Material cross-section 
  TMXsec *mxs = 0;
  // ((TMXsec*)((TGeoRCExtension*)mat−>
  //                        GetFWExtension()) −> GetUserObject() ) ;

  Int_t reactionId= -1;
  double kinEnergy= track.GetKineticEnergy()/CLHEP::GeV;
  // TEXsec* TMXsec::SampleInt(Int t part, Double_t en, Int_t &reac)
  TEXsec* elemXsec= mxs->SampleInt(fParticleId,
                                   kinEnergy,
                                   reactionId); // Output
  G4int indexElem = elemXsec->Index();
  
  TEFstate       **pElemFstate;   // Array of final state pointers per element
  // pElemFstate=   pTTabPhysMgr->GetElemFstate();
  pElemFstate=     theDataManager->GetElemFstate();
  
  Int_t nSecPart     = 0;  //number of secondary particles per reaction
  const Int_t *pid   = 0;  //GeantV particle codes [nSecPart]
  const Float_t *mom = 0;  //momentum vectors the secondaries [3*nSecPart]
  Float_t  energyFst = 0;  //energy at the fstate (Ekin of primary after the interc.)
  Float_t  kerma     = 0;  //released energy
  Float_t  weightFst = 0;  //weight of the fstate (just a dummy parameter now)
  Char_t   isSurv    = 0;  //is the primary survived the interaction
  
  isSurv = pElemFstate[indexElem]->SampleReac(fParticleId,
                                              reactionId,
                                              kinEnergy,
                                              nSecPart,
                                              weightFst,
                                              kerma,
                                              energyFst,
                                              pid,
                                              mom);
  
  G4int num = fSecParticles.size();
  if(num > 0) {
    fParticleChange->SetNumberOfSecondaries(num);
    G4double time = track.GetGlobalTime();
    G4double weight = fParticleChange->GetParentWeight();

    for(G4int is = 0 ; is < num ; ++is) {
      //create G4DaynamicParticle for EM processes
      const G4ParticleDefinition* particleDef= ParticleDefinition(fSecParticles[is]->id);

      G4ThreeVector secMomentum(fSecParticles[is]->px,
				fSecParticles[is]->py,
				fSecParticles[is]->px);

      G4DynamicParticle* dynamicParticle = 
         new G4DynamicParticle(particleDef,
                               secMomentum.unit(),
                               fSecParticles[is]->E);        

      G4Track* t = new G4Track(dynamicParticle, time, track.GetPosition());
      t->SetTouchableHandle(track.GetTouchableHandle());
      t->SetWeight(weight); 

      fParticleChange->AddSecondary(t);
    }
  }

  // Deal with change of the primary particle
  if( isSurv ) {
  
  }else{
     
  }

  
  return fParticleChange;
  */
}

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
