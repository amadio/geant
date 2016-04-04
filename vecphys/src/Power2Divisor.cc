#include "Power2Divisor.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST_DEVICE
Power2Divisor::Power2Divisor(int nmin, int nmax, int ndiv)
  : fNmin(nmin), fNmax(nmax), fNdiv(ndiv)
{
  fLowerBound = ldexp(1.,fNmin);
  fUpperBound = ldexp(1.,fNmax);
}

// following methods can be templated if necessary

VECCORE_CUDA_HOST_DEVICE
int Power2Divisor::GetNumberOfBins()
{
  //check for fNmax > fNmin
  return (fNmax-fNmin)*fNdiv+1;
}

VECCORE_CUDA_HOST_DEVICE
Precision Power2Divisor::GetBinPosition(int ibin)
{
  int exponent = Power2Exponent(ibin);
  //  int idiv = ibin & (fNdiv -1); // idiv=ibin%ndiv for any fNdiv = 2^n
  int idiv = ibin - fNdiv*Floor(ibin/fNdiv); //idiv=ibin%fNdiv (not restricted)
  return  ldexp(1.+ 1.*idiv/fNdiv,exponent);
}

VECCORE_CUDA_HOST_DEVICE
Precision Power2Divisor::GetBinSize(int ibin)
{
  int exponent = Power2Exponent(ibin);
  return ldexp(1.,exponent)/fNdiv;
}

VECCORE_CUDA_HOST_DEVICE
int Power2Divisor::Power2Exponent(int ibin) {
  return fNmin + ibin/fNdiv;
}

} // end namespace impl
} // end namespace vecphys
