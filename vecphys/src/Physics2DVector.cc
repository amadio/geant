#include "Physics2DVector.h"

// To-do : vectorize the Value method - any potential substitution for
//         FindBinLocationX/Y (binary search)?

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST_DEVICE Physics2DVector::Physics2DVector()
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

VECCORE_ATT_HOST_DEVICE
double Physics2DVector::Value(double x, double y)
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
  double x1 = xVector[idx];
  double x2 = xVector[idx+1];
  double y1 = yVector[idy];
  double y2 = yVector[idy+1];
  double v11= GetValue(idx,   idy);
  double v12= GetValue(idx+1, idy);
  double v21= GetValue(idx,   idy+1);
  double v22= GetValue(idx+1, idy+1);
  return ((y2 - y)*(v11*(x2 - x) + v12*(x - x1)) +
	  ((y - y1)*(v21*(x2 - x) + v22*(x - x1))))/((x2 - x1)*(y2 - y1));
}

VECCORE_ATT_HOST_DEVICE
void Physics2DVector::PutX(size_t idx, double val)
{
  xVector[idx] = val;
}

VECCORE_ATT_HOST_DEVICE
void Physics2DVector::PutY(size_t idy, double val)
{
  yVector[idy] = val;
}

VECCORE_ATT_HOST_DEVICE void
Physics2DVector::PutValue(size_t idx, size_t idy, double val)
{
  value[idy][idx] = val;
}

VECCORE_ATT_HOST_DEVICE
double Physics2DVector::GetValue(size_t idx, size_t idy)
{
  return value[idy][idx];
}

VECCORE_ATT_HOST_DEVICE
size_t Physics2DVector::FindBinLocationX(double z)
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

VECCORE_ATT_HOST_DEVICE
size_t Physics2DVector::FindBinLocationY(double z)
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

} // end namespace impl
} // end namespace vecphys
