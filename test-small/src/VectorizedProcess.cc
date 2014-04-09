//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "VectorizedProcess.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4GPILSelection.hh"
#include "G4SystemOfUnits.hh"

VectorizedProcess::VectorizedProcess(G4String processName, 
				     G4ProcessType processType) :
  G4WrapperProcess(processName,processType),
  particleChange(0) 
{
}

VectorizedProcess::~VectorizedProcess() {
}

G4VParticleChange* 
VectorizedProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // process PostStepDoIt for the original process
  particleChange = pRegProcess->PostStepDoIt(track, step);

  //  Print(step);

  return particleChange;
}

void VectorizedProcess::Print(const G4Step& step) 
{
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  std::cout << "VectorizedProcess ProcessName, Position, Energy, Nsecondary " 
	    << postStepPoint->GetProcessDefinedStep()->GetProcessName() << " " 
	    << postStepPoint->GetPosition() << " "  
	    << postStepPoint->GetKineticEnergy()/GeV  << " " 
	    << particleChange->GetNumberOfSecondaries() << std::endl;
}
