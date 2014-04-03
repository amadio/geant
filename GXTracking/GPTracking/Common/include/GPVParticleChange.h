#ifndef GPVParticleChange_h
#define GPVParticleChange_h 1

#include "GPConstants.h"
#include "GPStep.h"
#include "GPTrackStatus.h"
#include "GPThreeVector.h"
#include "GXTrack.h"

//class GPVParticleChange 

class GPVParticleChange 
{
public:
  FQUALIFIER GPVParticleChange();
  FQUALIFIER ~GPVParticleChange();

  FQUALIFIER void UpdateStepForAlongStep(GPStep* Step);
  FQUALIFIER void UpdateStepForPostStep(GPStep* Step);
  FQUALIFIER void UpdateStepInfo(GPStep* Step);

  FQUALIFIER void Initialize(GPStep* aStep);
  FQUALIFIER void InitializeTrueStepLength(GPStep* aStep);
  FQUALIFIER void InitializeLocalEnergyDeposit();
  FQUALIFIER void InitializeSecondaries();

  FQUALIFIER GPTrackStatus GetTrackStatus();
  FQUALIFIER void ProposeTrackStatus(GPTrackStatus aStatus);

  FQUALIFIER G4double GetLocalEnergyDeposit();
  FQUALIFIER void ProposeLocalEnergyDeposit(G4double anEnergyPart);

  FQUALIFIER G4double GetTrueStepLength();
  FQUALIFIER void  ProposeTrueStepLength(G4double truePathLength);

  FQUALIFIER G4int GetNumberOfSecondaries();
  FQUALIFIER GXTrack* GetSecondary(G4int anIndex);
  FQUALIFIER void AddSecondary(GXTrack& aSecondary);

  //G4ParticleChangeForLoss
  FQUALIFIER G4double GetProposedCharge();
  FQUALIFIER void SetProposedCharge(G4double theCharge);
  FQUALIFIER G4double GetProposedKineticEnergy();
  FQUALIFIER void SetProposedKineticEnergy(G4double proposedKinEnergy);
  FQUALIFIER GPThreeVector& GetProposedMomentumDirection();
  FQUALIFIER void SetProposedMomentumDirection(GPThreeVector& dir);
  FQUALIFIER void SetLowEnergyLimit(G4double elimit);

  FQUALIFIER
  void ParticleChangeForLoss_InitializeForAlongStep(GXTrack* track);

  FQUALIFIER
  void ParticleChangeForLoss_InitializeForPostStep(GXTrack* track);

  FQUALIFIER
  void ParticleChangeForLoss_UpdateStepForAlongStep(GPStep* pStep);

  FQUALIFIER
  void ParticleChangeForLoss_UpdateStepForPostStep(GPStep* pStep);

  //G4ParticleChangeForMSC
  FQUALIFIER void SetProposedPosition(GPThreeVector& P);
  FQUALIFIER GPThreeVector& GetProposedPosition();
  FQUALIFIER void ParticleChangeForMSC_Initialize(GXTrack* track);

  FQUALIFIER
  void ParticleChangeForMSC_UpdateStepForAlongStep(GPStep* pStep);

  FQUALIFIER
  void ParticleChangeForMSC_UpdateStepForPostStep(GPStep* pStep);

  FQUALIFIER
  void ParticleChangeForGamma_InitializeForPostStep(GXTrack* track);

  //Utility
  FQUALIFIER G4double CalculateKineticEnergy(GXTrack *track);
  //  FQUALIFIER void ClearSecondaries();
  FQUALIFIER void Clear();
 
private:

  GXTrack theListOfSecondaries[maximumNumberOfSecondaries];

  G4int theNumberOfSecondaries;
  GPTrackStatus theStatusChange;
  G4double theLocalEnergyDeposit;
  G4double theTrueStepLength;

  //G4ParticleChangeForLoss
  G4double proposedKinEnergy;
  G4double lowEnergyLimit;
  G4double currentCharge;

  //G4ParticleChangeForLoss+G4ParticleChangeForMSC
  GPThreeVector theMomentumDirection;

  //G4ParticleChangeForMSC
  GPThreeVector thePosition;

};

#endif
