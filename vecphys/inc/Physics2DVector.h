//a temporary interface to load and read SeltzerBerger data
//(adopted from G4Physics2DVector) - need to vectorize this utility class

#ifndef Physics2DVector_H
#define Physics2DVector_H 1

#include "base/Global.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class Physics2DVector
{
public:

  VECCORE_CUDA_HOST_DEVICE
  Physics2DVector();

  VECCORE_CUDA_HOST_DEVICE
  ~Physics2DVector(){};

  VECCORE_CUDA_HOST_DEVICE
  double Value(double x, double y);

  VECCORE_CUDA_HOST_DEVICE
  void PutX(size_t idx, double val);

  VECCORE_CUDA_HOST_DEVICE
  void PutY(size_t idy, double val);

  VECCORE_CUDA_HOST_DEVICE
  void PutValue(size_t idx, size_t idy, double val);

  VECCORE_CUDA_HOST_DEVICE
  double GetValue(size_t idx, size_t idy);

  VECCORE_CUDA_HOST_DEVICE
  size_t FindBinLocationX(double x);

  VECCORE_CUDA_HOST_DEVICE
  size_t FindBinLocationY(double y);

private:
  double  xVector[numberOfXNodes];
  double  yVector[numberOfYNodes];
  double  value[numberOfYNodes][numberOfXNodes];

};

} // end namespace impl
} // end namespace vecphys

#endif
