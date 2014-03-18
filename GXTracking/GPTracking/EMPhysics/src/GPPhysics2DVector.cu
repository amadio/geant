#include "GPPhysics2DVector.h"

FQUALIFIER GPPhysics2DVector::GPPhysics2DVector()
{
  // initialisation
  lastX = lastY = lastValue = 0.0;
  lastBinX = lastBinY = 0;  

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

// --------------------------------------------------------------
FQUALIFIER
void GPPhysics2DVector::ComputeValue(G4double xx, G4double yy)
{
  if(xx != lastBinX) {
    if(xx <= xVector[0]) {
      lastX = xVector[0];
      lastBinX = 0;
    } else if(xx >= xVector[numberOfXNodes-1]) {
      lastX = xVector[numberOfXNodes-1];
      lastBinX = numberOfXNodes-2;
    } else {
      lastX = xx;
      FindBinLocationX(xx);
    }
  }
  if(yy != lastBinY) {
    if(yy <= yVector[0]) {
      lastY = yVector[0];
      lastBinY = 0;
    } else if(yy >= yVector[numberOfYNodes-1]) {
      lastY = yVector[numberOfYNodes-1];
      lastBinY = numberOfYNodes-2;
    } else {
      lastY = yy;
      FindBinLocationY(yy);
    }
  }
  size_t idx  = lastBinX;
  size_t idy  = lastBinY;

  G4double x1 = xVector[idx];
  G4double x2 = xVector[idx+1];
  G4double y1 = yVector[idy];
  G4double y2 = yVector[idy+1];
  G4double x  = lastX;
  G4double y  = lastY;
  G4double v11= GetValue(idx,   idy);
  G4double v12= GetValue(idx+1, idy);
  G4double v21= GetValue(idx,   idy+1);
  G4double v22= GetValue(idx+1, idy+1);

  lastValue = 
    ((y2 - y)*(v11*(x2 - x) + v12*(x - x1)) + 
     ((y - y1)*(v21*(x2 - x) + v22*(x - x1))))/((x2 - x1)*(y2 - y1)); 
}

FQUALIFIER 
G4double GPPhysics2DVector::Value(G4double x, G4double y)
{
  if(x != lastX || y != lastY) { ComputeValue(x, y); }
  return lastValue;
}

FQUALIFIER
void GPPhysics2DVector::PutX(size_t idx, G4double val)
{
  xVector[idx] = val;
}

FQUALIFIER
void GPPhysics2DVector::PutY(size_t idy, G4double val)
{
  yVector[idy] = val;
}

FQUALIFIER void 
GPPhysics2DVector::PutValue(size_t idx, size_t idy, G4double val)
{
  value[idy][idx] = val;
}

FQUALIFIER
G4double GPPhysics2DVector::GetValue(size_t idx, size_t idy) const
{
  return value[idy][idx];
}

FQUALIFIER 
void GPPhysics2DVector::FindBinLocationX(G4double x)
{
  if(x < xVector[lastBinX] || x >=  xVector[lastBinX]) {  
    size_t lowerBound = 0;
    size_t upperBound = numberOfXNodes - 2;

    while (lowerBound <= upperBound) {
      size_t midBin = (lowerBound + upperBound)/2;
      if( x < xVector[midBin] ) { upperBound = midBin-1; }
      else                      { lowerBound = midBin+1; }
    }
    lastBinX = upperBound;
  }
}

FQUALIFIER
void GPPhysics2DVector::FindBinLocationY(G4double y)
{
  if(y < yVector[lastBinY] || y >=  yVector[lastBinY]) {  
    size_t lowerBound = 0;
    size_t upperBound = numberOfYNodes - 2;

    while (lowerBound <= upperBound) {
      size_t midBin = (lowerBound + upperBound)/2;
      if( y < yVector[midBin] ) { upperBound = midBin-1; }
      else                      { lowerBound = midBin+1; }
    }
    lastBinY = upperBound;
  }
}
