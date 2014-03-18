#include "GPStep.h"

FQUALIFIER
GPStep::GPStep()
  :  fTotalEnergyDeposit(0.0),
     fStepLength(0.)
{
}

FQUALIFIER
GPStep::~GPStep()
{
}

// Get/Set functions 
FQUALIFIER 
GPStepPoint& GPStep::GetPreStepPoint()
{ 
  return fpPreStepPoint; 
}

FQUALIFIER 
void GPStep::SetPreStepPoint(GPStepPoint& value)
{  
  fpPreStepPoint = value; 
}

FQUALIFIER 
GPStepPoint& GPStep::GetPostStepPoint()
{ 
  return fpPostStepPoint; 
}

FQUALIFIER 
void GPStep::SetPostStepPoint(GPStepPoint& value)
{ 
  fpPostStepPoint = value; 
}

FQUALIFIER
GPThreeVector GPStep::GetDeltaPosition()
{ 
  return GPThreeVector_sub(fpPostStepPoint.GetPosition(),
			   fpPreStepPoint.GetPosition()); 
}

FQUALIFIER
GPThreeVector GPStep::GetDeltaMomentum()
{ 
  return GPThreeVector_sub(fpPostStepPoint.GetMomentum(),
			   fpPreStepPoint.GetMomentum()); 
}

FQUALIFIER
G4double GPStep::GetDeltaEnergy()
{ 
  return fpPostStepPoint.GetKineticEnergy()
    - fpPreStepPoint.GetKineticEnergy(); 
}

FQUALIFIER 
G4double GPStep::GetStepLength()
{ 
  return fStepLength; 
}

FQUALIFIER 
void GPStep::SetStepLength(G4double value)
{ 
  fStepLength = value; 
}

FQUALIFIER 
G4double GPStep::GetTotalEnergyDeposit()
{ 
  return fTotalEnergyDeposit; 
}

FQUALIFIER 
void GPStep::SetTotalEnergyDeposit(G4double value)
{ 
  fTotalEnergyDeposit = value;   
}

FQUALIFIER 
void GPStep::AddTotalEnergyDeposit(G4double value)
{
  fTotalEnergyDeposit += value;   
}

FQUALIFIER 
void GPStep::ResetTotalEnergyDeposit()
{ 
  fTotalEnergyDeposit = 0.; 
}

FQUALIFIER 
void GPStep::CopyPostToPreStepPoint( )
{
  //This method is called at the beggining of each step 
  fpPreStepPoint = fpPostStepPoint;
  fpPostStepPoint.SetStepStatus(fUndefined);
}

FQUALIFIER
void GPStep::SetTrack(GXTrack* value)
{ 
  fpTrack = value; 
}

FQUALIFIER
void GPStep::InitializeStep(GXTrack* track )
{
  // Initialize G4Step attributes
  fStepLength = 0.;
  fTotalEnergyDeposit = 0.;
  fpTrack = track;
  //   fpTrack->SetStepLength(0.);
  
  // Initialize G4StepPoint attributes.
  // To avoid the circular dependency between G4Track, G4Step
  // and G4StepPoint, G4Step has to manage the copy actions.
  
  GPThreeVector pos = GPThreeVector_create(track->x,track->y,track->z);
  GPThreeVector dir = GPThreeVector_unit(
		      GPThreeVector_create(track->px,track->py,track->pz));

  G4double p = sqrt(track->px*track->px + track->py*track->py 
		    + track->pz*track->pz);
  G4double mass = electron_mass_c2*track->q*track->q;
  G4double kineticEnergy = p*p/(sqrt(p*p + mass*mass) + mass);
  
  fpPreStepPoint.SetPosition(pos);
  fpPreStepPoint.SetMomentumDirection(dir);
  fpPreStepPoint.SetKineticEnergy(kineticEnergy);
  fpPreStepPoint.SetSafety(0);
  fpPreStepPoint.SetStepStatus(fUndefined);
  fpPreStepPoint.SetMass(mass);    
  fpPreStepPoint.SetCharge(track->q); 
  
  fpPostStepPoint = fpPreStepPoint;
}

FQUALIFIER
void GPStep::UpdateTrack()
{ 
  //  fpTrack->SetPosition(fpPostStepPoint->GetPosition());
  fpTrack->x = (fpPostStepPoint.GetPosition()).x;
  fpTrack->y = (fpPostStepPoint.GetPosition()).y;
  fpTrack->z = (fpPostStepPoint.GetPosition()).z;

  fpTrack->px = (fpPostStepPoint.GetMomentum()).x;
  fpTrack->py = (fpPostStepPoint.GetMomentum()).y;
  fpTrack->pz = (fpPostStepPoint.GetMomentum()).z;

  fpTrack->E = fpPostStepPoint.GetKineticEnergy();
  fpTrack->q = fpPostStepPoint.GetCharge();
  fpTrack->s = fStepLength;

  //  fpTrack->SetNextTouchableHandle(fpPostStepPoint->GetTouchableHandle());
  //  fpTrack->SetWeight(fpPostStepPoint->GetWeight());
  
}
