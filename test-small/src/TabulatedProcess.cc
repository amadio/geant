//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "TabulatedProcess.hh"
#include "TabulatedDataManager.hh"
#include "TabulatedPhysicsNameSpace.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4GPILSelection.hh"
#include "G4SystemOfUnits.hh"

/*
#include "TSystem.h"
#include <TPartIndex.h>
#include <TEXsec.h>
#include <TEFstate.h>
*/

TabulatedProcess::TabulatedProcess(G4String processName, 
				   G4ProcessType processType) :
  G4WrapperProcess(processName,processType),
  particleChange(0)
{
  theDataManager = TabulatedDataManager::Instance();
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
  return pRegProcess->PostStepGetPhysicalInteractionLength(track,
						    previousStepSize,
						    condition);

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

G4double TabulatedProcess::MeanFreePath(const G4Track& track)
{
  G4double energy = track.GetKineticEnergy()/GeV;

  //material index and the number of particles with reactions
  G4int imat = 1;  //temporary
  G4int npart = 1; //temporary

  G4double preStepLambda = 
    theDataManager->GetInteractionLength(imat,npart,energy);
  return preStepLambda;
}

G4VParticleChange* 
TabulatedProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // process PostStepDoIt for the original process
  particleChange = pRegProcess->PostStepDoIt(track, step);

  //implement

  return particleChange;
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
	    << particleChange->GetNumberOfSecondaries() << std::endl;
}

