#ifndef GXPhysics2DVector_H
#define GXPhysics2DVector_H 1

#include "GPTypeDef.h"
#include "GPConstants.h"

class GXPhysics2DVector 
{
public: 

  FQUALIFIER GXPhysics2DVector();
  FQUALIFIER ~GXPhysics2DVector(){};

  FQUALIFIER G4double Value(G4double x, G4double y);

  FQUALIFIER void PutX(size_t idx, G4double val);
  FQUALIFIER void PutY(size_t idy, G4double val);
  FQUALIFIER void PutValue(size_t idx, size_t idy, G4double val);
  FQUALIFIER G4double GetValue(size_t idx, size_t idy);

  FQUALIFIER size_t FindBinLocationX(G4double x);
  FQUALIFIER size_t FindBinLocationY(G4double y);

private:
  G4double  xVector[numberOfXNodes];
  G4double  yVector[numberOfYNodes];
  G4double  value[numberOfYNodes][numberOfXNodes];

};

#endif
