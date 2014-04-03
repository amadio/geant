#ifndef GPPhysics2DVector_H
#define GPPhysics2DVector_H 1

#include "GPTypeDef.h"
#include "GPConstants.h"

class GPPhysics2DVector 
{
public: 

  FQUALIFIER GPPhysics2DVector();
  FQUALIFIER ~GPPhysics2DVector(){};

  FQUALIFIER void ComputeValue(G4double x, G4double y);
  FQUALIFIER G4double Value(G4double x, G4double y);

  FQUALIFIER void PutX(size_t idx, G4double val);
  FQUALIFIER void PutY(size_t idy, G4double val);
  FQUALIFIER void PutValue(size_t idx, size_t idy, G4double value);
  FQUALIFIER G4double GetValue(size_t idx, size_t idy) const;

  FQUALIFIER void FindBinLocationX(G4double x);
  FQUALIFIER void FindBinLocationY(G4double y);

private:
  G4double  xVector[numberOfXNodes];
  G4double  yVector[numberOfYNodes];
  G4double  value[numberOfYNodes][numberOfXNodes];

  //cache;
  G4double lastX;
  G4double lastY;
  G4double lastValue;
  size_t   lastBinX;
  size_t   lastBinY;
};

#endif
