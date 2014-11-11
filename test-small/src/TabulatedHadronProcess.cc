//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "TabulatedHadronProcess.hh"
#include "TabulatedDataManager.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4GPILSelection.hh"
#include "G4SystemOfUnits.hh"

TabulatedHadronProcess::TabulatedHadronProcess(G4String processName, 
				   G4ProcessType processType) :
  G4WrapperProcess(processName,processType),
  particleChange(0),
  theDataManager(TabulatedDataManager::Instance())
{
}

TabulatedHadronProcess::~TabulatedHadronProcess() {
}

G4double
TabulatedHadronProcess::PostStepGetPhysicalInteractionLength(const G4Track& track,
						   G4double previousStepSize,
						   G4ForceCondition* condition)
{
  // processPostStepGetPhysicalInteractionLength for the original process
  return 
    pRegProcess->PostStepGetPhysicalInteractionLength(track, previousStepSize,
						      condition);
  //replace by the VP tabulated physics

}

G4VParticleChange* 
TabulatedHadronProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // process PostStepDoIt for the original process
  particleChange = pRegProcess->PostStepDoIt(track, step);

  //replace by the VP tabulated physics
  //  Print(step);

  return particleChange;
}

void TabulatedHadronProcess::Print(const G4Step& step) 
{
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  std::cout << "TabulatedHadronProcess ProcName, E, dE, x, StepLength Nsec " 
	    << postStepPoint->GetProcessDefinedStep()->GetProcessName() << " " 
	    << postStepPoint->GetKineticEnergy()/MeV  << " " 
	    << step.GetTotalEnergyDeposit()/MeV  << " " 
	    << postStepPoint->GetPosition() << " "  
	    << step.GetStepLength() << " " 
	    << particleChange->GetNumberOfSecondaries() << std::endl;
}
