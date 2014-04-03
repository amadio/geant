#ifndef GPStepPoint_H
#define GPStepPoint_H 1

#include "GPThreeVector.h"
#include "GPStepStatus.h" 

class GPStepPoint
{
public:

  FQUALIFIER GPStepPoint();
  FQUALIFIER ~GPStepPoint(){}

  FQUALIFIER GPStepPoint(const GPStepPoint& right);
  FQUALIFIER GPStepPoint& operator=(const GPStepPoint &);   

  FQUALIFIER GPThreeVector& GetPosition() ;
  FQUALIFIER void SetPosition(GPThreeVector& aValue);
  FQUALIFIER void AddPosition(GPThreeVector& aValue);

  FQUALIFIER GPThreeVector& GetMomentumDirection() ;
  FQUALIFIER void SetMomentumDirection(GPThreeVector& aValue);
  FQUALIFIER void AddMomentumDirection(GPThreeVector& aValue);
    
  FQUALIFIER GPThreeVector GetMomentum();
  FQUALIFIER G4double GetTotalEnergy();

  FQUALIFIER G4double GetKineticEnergy();
  FQUALIFIER void SetKineticEnergy(G4double aValue);
  FQUALIFIER void AddKineticEnergy(G4double aValue);
  
  FQUALIFIER G4double GetSafety();
  FQUALIFIER void SetSafety(G4double aValue);
  
  FQUALIFIER GPStepStatus GetStepStatus();
  FQUALIFIER void SetStepStatus(GPStepStatus aValue);
  
  FQUALIFIER G4double GetMass();
  FQUALIFIER void SetMass(G4double value);
  
  FQUALIFIER G4double GetCharge();
  FQUALIFIER void SetCharge(G4double value);
  
private:
  GPThreeVector fPosition;
  GPThreeVector fMomentumDirection;
  G4double fKineticEnergy;
  G4double fMass;
  G4double fCharge;
  G4double fSafety;
  GPStepStatus fStepStatus;

};

#endif
