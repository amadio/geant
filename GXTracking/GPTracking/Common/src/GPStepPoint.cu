#include "GPStepPoint.h"

FQUALIFIER
GPStepPoint::GPStepPoint() :
 fKineticEnergy(0.),fMass(0.), fCharge(0.),fStepStatus(fUndefined)
{
}

FQUALIFIER
GPStepPoint::GPStepPoint(const GPStepPoint &right) :
  fPosition(right.fPosition),
  fMomentumDirection(right.fMomentumDirection),
  fKineticEnergy(right.fKineticEnergy),
  fMass(right.fMass),
  fCharge(right.fCharge),
  fStepStatus(right.fStepStatus)
{
}

FQUALIFIER
GPStepPoint& GPStepPoint::operator=(const GPStepPoint &right)
{
  if (this != &right) {
    fPosition     = right.fPosition;
    fMomentumDirection = right.fMomentumDirection;
    fKineticEnergy = right.fKineticEnergy;
    fMass         = right.fMass;
    fCharge       = right.fCharge;
    fStepStatus   = right.fStepStatus;
  }
  return *this;
}

FQUALIFIER 
GPThreeVector& GPStepPoint::GetPosition()
{ 
  return fPosition; 
}

FQUALIFIER 
void GPStepPoint::SetPosition(GPThreeVector& aValue)
{ 
  fPosition = aValue; 
}

FQUALIFIER 
void GPStepPoint::AddPosition(GPThreeVector& aValue)
{ 
  fPosition = GPThreeVector_add(fPosition,aValue); 
} 
    
FQUALIFIER 
GPThreeVector& GPStepPoint::GetMomentumDirection()
{ 
  return fMomentumDirection; 
}

FQUALIFIER 
void GPStepPoint::SetMomentumDirection(GPThreeVector& aValue)
{ 
  fMomentumDirection = aValue; 
}

FQUALIFIER
void GPStepPoint::AddMomentumDirection(GPThreeVector& aValue)
{ 
  fMomentumDirection = GPThreeVector_add(fMomentumDirection,aValue);
} 
    
FQUALIFIER 
GPThreeVector GPStepPoint::GetMomentum()
{ 
  G4double tMomentum = sqrt(fKineticEnergy*fKineticEnergy +
			    2*fKineticEnergy*fMass);

  return GPThreeVector_create(fMomentumDirection.x*tMomentum,
			      fMomentumDirection.y*tMomentum,
			      fMomentumDirection.z*tMomentum);
}

FQUALIFIER 
G4double GPStepPoint::GetTotalEnergy()
{ 
  return fKineticEnergy + fMass; 
} 

FQUALIFIER
G4double GPStepPoint::GetKineticEnergy()
{ 
  return fKineticEnergy; 
}

FQUALIFIER
void GPStepPoint::SetKineticEnergy(G4double aValue)
{ 
  fKineticEnergy = aValue; 
}

FQUALIFIER
void GPStepPoint::AddKineticEnergy(G4double aValue)
{ 
  fKineticEnergy += aValue; 
}

FQUALIFIER
GPStepStatus GPStepPoint::GetStepStatus()
{ 
  return fStepStatus; 
}

FQUALIFIER
void GPStepPoint::SetStepStatus(GPStepStatus aValue)
{ 
  fStepStatus = aValue; 
}

FQUALIFIER
G4double GPStepPoint::GetMass()
{ 
  return fMass; 
}

FQUALIFIER
void GPStepPoint::SetMass(G4double value)
{ 
  fMass = value; 
}

FQUALIFIER
G4double GPStepPoint::GetCharge()
{ 
  return fCharge; 
}

FQUALIFIER
void GPStepPoint::SetCharge(G4double value)
{ 
  fCharge = value;
}

FQUALIFIER
G4double GPStepPoint::GetSafety()
{ 
  return fSafety; 
}

FQUALIFIER
void GPStepPoint::SetSafety(G4double aValue)
{ 
  fSafety = aValue; 
}
