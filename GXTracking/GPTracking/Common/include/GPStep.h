#ifndef G4Step_h
#define G4Step_h 1

#include "GPThreeVector.h"  
#include "GPStepPoint.h"
#include "GXTrack.h"

class GPStep
{

public:

  FQUALIFIER GPStep();
  FQUALIFIER ~GPStep();

  FQUALIFIER GPStep(const GPStep& );
  FQUALIFIER GPStep & operator=(const GPStep &);   

  // step points 
  FQUALIFIER GPStepPoint& GetPreStepPoint();
  FQUALIFIER void SetPreStepPoint(GPStepPoint& value);

  FQUALIFIER GPStepPoint& GetPostStepPoint();
  FQUALIFIER void SetPostStepPoint(GPStepPoint& value);

  FQUALIFIER GPThreeVector GetDeltaPosition();
  FQUALIFIER GPThreeVector GetDeltaMomentum();
  FQUALIFIER G4double GetDeltaEnergy();

  FQUALIFIER G4double GetStepLength();
  FQUALIFIER void SetStepLength(G4double value);

  FQUALIFIER G4double GetTotalEnergyDeposit();
  FQUALIFIER void SetTotalEnergyDeposit(G4double value);

  FQUALIFIER void AddTotalEnergyDeposit(G4double value);
  FQUALIFIER void ResetTotalEnergyDeposit();

  FQUALIFIER void SetTrack(GXTrack* value);
  FQUALIFIER void InitializeStep( GXTrack* aValue );
  FQUALIFIER void UpdateTrack();
  FQUALIFIER void CopyPostToPreStepPoint( );
  
private:

  GPStepPoint fpPreStepPoint;
  GPStepPoint fpPostStepPoint;
  GXTrack *fpTrack;
  G4double fStepLength;
  G4double fTotalEnergyDeposit;
};

#endif
