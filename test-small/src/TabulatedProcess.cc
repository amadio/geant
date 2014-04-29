//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "TabulatedProcess.hh"
#include "TabulatedDataManager.hh"
#include "TabulatedPhysicsNameSpace.hh"

#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
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
#include <TEXsec.h>
#include <TEFstate.h>
*/

TabulatedProcess::TabulatedProcess(G4String processName, 
				   G4ProcessType processType) :
  G4WrapperProcess(processName,processType),
  fMaterialIndex(-1),
  fReaction(kTotal),
  fParticleChange(0)
{
  fSecDefinition = G4Geantino::Geantino();
  theDataManager = 0;
  fSecParticles.reserve(5);
}

TabulatedProcess::TabulatedProcess(G4String processName, 
				   G4ProcessType processType,
				   G5proc reactionIndex) :
  G4WrapperProcess(processName,processType),
  fMaterialIndex(-1),
  fReaction(reactionIndex),
  fParticleChange(0)
{
  fSecDefinition = G4Geantino::Geantino();
  theDataManager = TabulatedDataManager::Instance();
  fSecParticles.reserve(5);
}

TabulatedProcess::~TabulatedProcess() 
{
}

G4double
TabulatedProcess::PostStepGetPhysicalInteractionLength(const G4Track& track,
						   G4double previousStepSize,
						   G4ForceCondition* condition)
{
  // Call to ensure that the process is correctly initialised - temporary
  pRegProcess->PostStepGetPhysicalInteractionLength(track,
						    previousStepSize,
						    condition);

  SetupForMaterial(track);

  if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
    // beggining of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else if ( previousStepSize > 0.0) {
    // subtract NumberOfInteractionLengthLeft 
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
  } else {
    // zero step
    //  DO NOTHING
  }

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean free path
  currentInteractionLength = MeanFreePath(track);

  G4double value;
  if (currentInteractionLength <DBL_MAX) {
    value = theNumberOfInteractionLengthLeft * currentInteractionLength;
  } else {
    value = DBL_MAX;
  }
  return value;
}

G4double TabulatedProcess::SetupForMaterial(const G4Track& track)
{
  //get the material index for the current volume
  fMaterialIndex = 1;  //temporary 
}

G4double TabulatedProcess::MeanFreePath(const G4Track& track)
{
  G4double preStepLambda = 
    theDataManager->GetInteractionLength(fMaterialIndex,track);

  return preStepLambda;
}

G4VParticleChange* 
TabulatedProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // process PostStepDoIt for the original process
  fParticleChange = pRegProcess->PostStepDoIt(track, step);

  //implementation for the wrapper process

  theNumberOfInteractionLengthLeft = -1.0;

  fSecParticles.clear();
  theDataManager->SampleSecondaries(&fSecParticles,fMaterialIndex,
				    &track,fReaction);

  G4int num = fSecParticles.size();
  if(num > 0) {
    fParticleChange->SetNumberOfSecondaries(num);
    G4double time = track.GetGlobalTime();
    G4double weight = fParticleChange->GetParentWeight();

    for(G4int is = 0 ; is < num ; ++is) {
      //create G4DaynamicParticle for EM processes
      ResetSecondaryDefinition(fSecParticles[is]->id);
      G4ThreeVector secMomentum(fSecParticles[is]->px,
				fSecParticles[is]->py,
				fSecParticles[is]->px);

      G4DynamicParticle* dynamicParticle = 
	new G4DynamicParticle(fSecDefinition,secMomentum.unit(),
			      fSecParticles[is]->E);        

      G4Track* t = new G4Track(dynamicParticle, time, track.GetPosition());
      t->SetTouchableHandle(track.GetTouchableHandle());
      t->SetWeight(weight); 

      fParticleChange->AddSecondary(t);
    }
  }

  //changes for the primary particle

  return fParticleChange;
}

void TabulatedProcess::ResetSecondaryDefinition(G4int pid) 
{
  //only for gamma/electron/photon for the standard EM processes

  switch (pid) {
  case 22:
    fSecDefinition = G4Gamma::Gamma(); 
    break;
  case 11:
    fSecDefinition = G4Electron::Electron(); 
    break;
  case -11:
    fSecDefinition = G4Positron::Positron(); 
    break;
  default:
    break;
  }
}

void TabulatedProcess::Print(const G4Step& step) 
{
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  std::cout << "TabulatedProcess ProcName, E, dE, x, StepLength Nsec " 
	    << postStepPoint->GetProcessDefinedStep()->GetProcessName() << " " 
	    << postStepPoint->GetKineticEnergy()/MeV  << " " 
	    << step.GetTotalEnergyDeposit()/MeV  << " " 
	    << postStepPoint->GetPosition() << " "  
	    << step.GetStepLength() << " " 
	    << fParticleChange->GetNumberOfSecondaries() << std::endl;
}

