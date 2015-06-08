//a temporary interface to load and read SeltzerBerger data 
//(adopted from G4Physics2DVector) - need to vectorize this utility class

#ifndef Physics2DVector_H
#define Physics2DVector_H 1

#include "core/base/Global.h"
#include "GUConstants.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

class Physics2DVector  
{
public: 

  VECPHYS_FUNC_QUALIFIER 
  Physics2DVector();

  VECPHYS_FUNC_QUALIFIER 
  ~Physics2DVector(){};

  VECPHYS_FUNC_QUALIFIER 
  double Value(double x, double y);

  VECPHYS_FUNC_QUALIFIER 
  void PutX(size_t idx, double val);

  VECPHYS_FUNC_QUALIFIER 
  void PutY(size_t idy, double val);

  VECPHYS_FUNC_QUALIFIER 
  void PutValue(size_t idx, size_t idy, double val);

  VECPHYS_FUNC_QUALIFIER 
  double GetValue(size_t idx, size_t idy);

  VECPHYS_FUNC_QUALIFIER 
  size_t FindBinLocationX(double x);

  VECPHYS_FUNC_QUALIFIER 
  size_t FindBinLocationY(double y);

private:
  double  xVector[numberOfXNodes];
  double  yVector[numberOfYNodes];
  double  value[numberOfYNodes][numberOfXNodes];

};

} // end namespace impl
} // end namespace vecphys

#endif
