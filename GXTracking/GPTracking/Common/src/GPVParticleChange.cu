#include "GPVParticleChange.h"
#include "GPSystemOfUnits.h"
#include <assert.h>

FQUALIFIER
GPVParticleChange::GPVParticleChange()
  :theNumberOfSecondaries(0),
   theStatusChange(fAlive),
   theLocalEnergyDeposit(0.0),
   theTrueStepLength(0.0)
{
  //G4ParticleChangeForLoss
  proposedKinEnergy = 0.;
  lowEnergyLimit = 1.0*eV; 
  currentCharge = 0.;
}

FQUALIFIER
GPVParticleChange::~GPVParticleChange() {
}

FQUALIFIER
void GPVParticleChange::AddSecondary(GXTrack &aTrack)
{
  // add a secondary after size check
  if (theNumberOfSecondaries < maximumNumberOfSecondaries) {
    theListOfSecondaries[theNumberOfSecondaries] = aTrack;
    theNumberOfSecondaries++;
  } else {
    assert(0);
  }
}
 
// Virtual methods for updating GPStep 
FQUALIFIER
void GPVParticleChange::UpdateStepInfo(GPStep* pStep)
{
  // Update the GPStep specific attributes
  pStep->SetStepLength( theTrueStepLength );
  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );
}

FQUALIFIER
void GPVParticleChange::UpdateStepForAlongStep(GPStep* Step)
{
  UpdateStepInfo(Step);
}

FQUALIFIER
void GPVParticleChange::UpdateStepForPostStep(GPStep* Step)
{
  UpdateStepInfo(Step);
}

//G4VParticleChange.icc : Set/Get FQUALIFIER functions

FQUALIFIER
GXTrack* GPVParticleChange::GetSecondary(G4int anIndex)
{
  return &(theListOfSecondaries[anIndex]);
}

FQUALIFIER 
G4int GPVParticleChange::GetNumberOfSecondaries()
{
  return theNumberOfSecondaries;
}

FQUALIFIER 
GPTrackStatus GPVParticleChange::GetTrackStatus()
{
  return theStatusChange;
}

FQUALIFIER  
void GPVParticleChange::ProposeTrackStatus(GPTrackStatus aStatus)
{
  theStatusChange = aStatus;
}

FQUALIFIER 
G4double GPVParticleChange::GetLocalEnergyDeposit()
{
  return theLocalEnergyDeposit;
}

FQUALIFIER 
void GPVParticleChange::ProposeLocalEnergyDeposit(G4double anEnergyPart)
{
  theLocalEnergyDeposit = anEnergyPart;
}

FQUALIFIER 
G4double GPVParticleChange::GetTrueStepLength()
{
  return theTrueStepLength;
}

FQUALIFIER 
void GPVParticleChange::ProposeTrueStepLength(G4double aLength)
{
  theTrueStepLength = aLength;
}

FQUALIFIER 
void GPVParticleChange::InitializeLocalEnergyDeposit()
{  
  // clear theLocalEnergyDeposited   
  theLocalEnergyDeposit = 0.0;
}

FQUALIFIER 
void GPVParticleChange::InitializeTrueStepLength(GPStep* aStep)
{
  theTrueStepLength = aStep->GetStepLength();
}

FQUALIFIER 
void GPVParticleChange::InitializeSecondaries()
{
  theNumberOfSecondaries = 0;
}

FQUALIFIER 
void GPVParticleChange::Initialize(GPStep* aStep)
{
  InitializeLocalEnergyDeposit();
  InitializeTrueStepLength(aStep);
  InitializeSecondaries();
}

// GPVParticleChange

FQUALIFIER 
G4double GPVParticleChange::GetProposedKineticEnergy()
{
  return proposedKinEnergy;
}

FQUALIFIER 
void GPVParticleChange::SetProposedKineticEnergy(G4double energy)
{
  proposedKinEnergy = energy;
}

FQUALIFIER 
G4double GPVParticleChange::GetProposedCharge()
{
  return currentCharge;
}

FQUALIFIER 
void GPVParticleChange::SetProposedCharge(G4double theCharge)
{
  currentCharge = theCharge;
}

FQUALIFIER
GPThreeVector& GPVParticleChange::GetProposedMomentumDirection() 
{
  return theMomentumDirection;
}

FQUALIFIER
void GPVParticleChange::SetProposedMomentumDirection(GPThreeVector& dir)
{
  theMomentumDirection = dir;
}

FQUALIFIER
void GPVParticleChange::SetProposedPosition(GPThreeVector& P)
{
  thePosition = P;
}

FQUALIFIER
GPThreeVector& GPVParticleChange::GetProposedPosition()
{
  return thePosition;
}

FQUALIFIER 
void GPVParticleChange::SetLowEnergyLimit(G4double elimit)
{
  lowEnergyLimit = elimit;
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForLoss_InitializeForAlongStep(
						         GXTrack* track)
{
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries();
  proposedKinEnergy = CalculateKineticEnergy(track);
  //  currentCharge = track.GetDynamicParticle()->GetCharge();
  currentCharge = track->q;
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForLoss_InitializeForPostStep(
                                                       GXTrack* track)
{
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries();
  proposedKinEnergy = CalculateKineticEnergy(track);
  //  currentCharge = track.GetDynamicParticle()->GetCharge();
  currentCharge = track->q;
  //  theMomentumDirection = track.GetMomentumDirection();
  theMomentumDirection = GPThreeVector_unit(
		       GPThreeVector_create(track->px,track->py,track->pz));
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForLoss_UpdateStepForAlongStep(
                                                       GPStep* pStep)
{
  GPStepPoint pPostStepPoint = pStep->GetPostStepPoint();

  // accumulate change of the kinetic energy
  G4double kinEnergy = pPostStepPoint.GetKineticEnergy() +
    (proposedKinEnergy - (pStep->GetPreStepPoint()).GetKineticEnergy());

  // update kinetic energy and charge
  if (kinEnergy < lowEnergyLimit) {
    theLocalEnergyDeposit += kinEnergy;
    kinEnergy = 0.0;
  } else {
    pPostStepPoint.SetCharge( currentCharge );
  }
  pPostStepPoint.SetKineticEnergy( kinEnergy );

  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );

  //@@@also put pPostStepPoint back to pStep
  pStep->SetPostStepPoint(pPostStepPoint);

  //  return pStep;
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForLoss_UpdateStepForPostStep(
                                                       GPStep* pStep)
{
  GPStepPoint pPostStepPoint = pStep->GetPostStepPoint();
  pPostStepPoint.SetCharge( currentCharge );
  pPostStepPoint.SetMomentumDirection( theMomentumDirection );
  pPostStepPoint.SetKineticEnergy( proposedKinEnergy );

  pStep->AddTotalEnergyDeposit( theLocalEnergyDeposit );

  //@@@also put pPostStepPoint back to pStep
  pStep->SetPostStepPoint(pPostStepPoint);

  //  return pStep;
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForMSC_Initialize(GXTrack* track)
{
  theMomentumDirection = GPThreeVector_unit(
		       GPThreeVector_create(track->px,track->py,track->pz));
  thePosition = GPThreeVector_create(track->x,track->y,track->z);
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForMSC_UpdateStepForAlongStep(
							     GPStep* pStep)
{
  //  Update the G4Step specific attributes
  pStep->SetStepLength(theTrueStepLength);
  //  theStatusChange = pStep->GetTrack()->GetTrackStatus();

  // Multiple scattering calculates the final state of the particle
  GPStepPoint pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint.SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint.SetPosition( thePosition );

  //@@@put pPostStepPoint back to pStep
  pStep->SetPostStepPoint(pPostStepPoint);
  //  return pStep;
}

FQUALIFIER
void GPVParticleChange::ParticleChangeForMSC_UpdateStepForPostStep(
                                                              GPStep* pStep)
{
  // Multiple scattering calculates the final state of the particle
  GPStepPoint pPostStepPoint = pStep->GetPostStepPoint();

  // update  momentum direction
  pPostStepPoint.SetMomentumDirection(theMomentumDirection);

  // update position
  pPostStepPoint.SetPosition( thePosition );

  //@@@put pPostStepPoint back to pStep
  pStep->SetPostStepPoint(pPostStepPoint);

  //  Update the G4Step specific attributes
  //  return pStep;
}

FQUALIFIER void 
GPVParticleChange::ParticleChangeForGamma_InitializeForPostStep(GXTrack* track)
{
  theLocalEnergyDeposit = 0.0;
  InitializeSecondaries();
  proposedKinEnergy = CalculateKineticEnergy(track);
  theMomentumDirection = GPThreeVector_unit(
		       GPThreeVector_create(track->px,track->py,track->pz));
}


// Utility functions

FQUALIFIER
G4double GPVParticleChange::CalculateKineticEnergy(GXTrack *track) {
  G4double p = sqrt(track->px*track->px + track->py*track->py 
		  + track->pz*track->pz);
  G4double mass = electron_mass_c2*track->q*track->q;
  G4double ekin = p*p/(sqrt(p*p + mass*mass) + mass);
  // the above expression is numerically more stable for both small and 
  // large momenta - syjun
  //  G4double ekin = sqrt(p*p + mass*mass) - mass;  // this is simpler!
  return ekin;
}

FQUALIFIER
void GPVParticleChange::Clear() {
  theNumberOfSecondaries = 0;
}
