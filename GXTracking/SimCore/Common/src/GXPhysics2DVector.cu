#include "GXPhysics2DVector.h"

FQUALIFIER GXPhysics2DVector::GXPhysics2DVector()
{
  for(size_t j = 0; j<numberOfYNodes; ++j) {
    yVector[j] = 0.0;
  }

  for(size_t i=0; i<numberOfXNodes; ++i) {
    xVector[i] = 0.0;
    for(size_t j=0; j<numberOfYNodes; ++j) {
      value[j][i] = 0.0;
    }
  }
};

FQUALIFIER 
G4double GXPhysics2DVector::Value(G4double x, G4double y)
{
  // no interpolation outside the table
  if(x < xVector[0]) { 
    x = xVector[0]; 
  } else if(x > xVector[numberOfXNodes - 1]) { 
    x = xVector[numberOfXNodes - 1]; 
  }
  if(y < yVector[0]) { 
    y = yVector[0]; 
  } else if(y > yVector[numberOfYNodes - 1]) { 
    y = yVector[numberOfYNodes - 1]; 
  }

  // find bins
  size_t idx = FindBinLocationX(x);
  size_t idy = FindBinLocationY(y);

  // interpolate
  //  if(useBicubic) {
  //    return BicubicInterpolation(x, y, idx, idy);
  //  } else {
  // useBicubic = false
  G4double x1 = xVector[idx];
  G4double x2 = xVector[idx+1];
  G4double y1 = yVector[idy];
  G4double y2 = yVector[idy+1];
  G4double v11= GetValue(idx,   idy);
  G4double v12= GetValue(idx+1, idy);
  G4double v21= GetValue(idx,   idy+1);
  G4double v22= GetValue(idx+1, idy+1);
  return ((y2 - y)*(v11*(x2 - x) + v12*(x - x1)) + 
	  ((y - y1)*(v21*(x2 - x) + v22*(x - x1))))/((x2 - x1)*(y2 - y1)); 
  //  }
}

FQUALIFIER
void GXPhysics2DVector::PutX(size_t idx, G4double val)
{
  xVector[idx] = val;
}

FQUALIFIER
void GXPhysics2DVector::PutY(size_t idy, G4double val)
{
  yVector[idy] = val;
}

FQUALIFIER void 
GXPhysics2DVector::PutValue(size_t idx, size_t idy, G4double val)
{
  value[idy][idx] = val;
}

FQUALIFIER
G4double GXPhysics2DVector::GetValue(size_t idx, size_t idy)
{
  return value[idy][idx];
}

FQUALIFIER 
size_t GXPhysics2DVector::FindBinLocationX(G4double z)
{
  size_t id = 0;
  if(z < xVector[1]) { 
    id = 0; 
  } 
  else if(z >= xVector[numberOfXNodes-2]) { 
    id = numberOfXNodes - 2; 
  } 
  else {
    size_t lowerBound = 0;
    size_t upperBound = numberOfXNodes - 2;

    while (lowerBound <= upperBound) {
      size_t midBin = (lowerBound + upperBound)/2;
      if( z < xVector[midBin] ) { upperBound = midBin-1; }
      else                      { lowerBound = midBin+1; }
    }
    id = upperBound;
  }
  return id;
}

FQUALIFIER 
size_t GXPhysics2DVector::FindBinLocationY(G4double z)
{
  size_t id = 0;
  if(z < yVector[1]) { 
    id = 0; 
  } 
  else if(z >= yVector[numberOfYNodes-2]) { 
    id = numberOfYNodes - 2; 
  } 
  else {
    size_t lowerBound = 0;
    size_t upperBound = numberOfYNodes - 2;

    while (lowerBound <= upperBound) {
      size_t midBin = (lowerBound + upperBound)/2;
      if( z < yVector[midBin] ) { upperBound = midBin-1; }
      else                      { lowerBound = midBin+1; }
    }
    id = upperBound;
  }
  return id;
}
