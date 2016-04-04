#ifndef POWER2DIVISOR_H
#define POWER2DIVISOR_H 1

#include "base/Global.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

class Power2Divisor
{
public:
  VECCORE_CUDA_HOST_DEVICE
  Power2Divisor(int nmin, int nmax, int ndiv);

  VECCORE_CUDA_HOST_DEVICE
  int GetNumberOfBins();

  VECCORE_CUDA_HOST_DEVICE
  Real_t GetLowerBound() { return fLowerBound; }

  VECCORE_CUDA_HOST_DEVICE
  Real_t GetUpperBound() { return fUpperBound; }

  VECCORE_CUDA_HOST_DEVICE
  Real_t GetBinPosition(int ibin);

  VECCORE_CUDA_HOST_DEVICE
  Real_t GetBinSize(int ibin);

  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void GetBinAndFraction(typename Backend::Double_v x,
                         typename Backend::Index_t& ibin,
                         typename Backend::Double_v& frac);

private:
  VECCORE_CUDA_HOST_DEVICE
  int Power2Exponent(int ibin);

private:
  int fNmin;              // minimum of the power2 exponent
  int fNmax;              // maximum of the power2 exponent
  int fNdiv;              // number of equal divisions between 2^e and 2^{e+1}

  Real_t fLowerBound;   // lower bound of the power2 range
  Real_t fUpperBound;   // upper bound of the power2 range
};

template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
void Power2Divisor::GetBinAndFraction(typename Backend::Double_v x,
                                      typename Backend::Index_t& ibin,
                                      typename Backend::Double_v& frac)
{
  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Index_t Index_t;
  typedef typename Backend::Double_v Double_v;

  Int_t    exponent;
  Double_v mantissa = frexp (x, &exponent); // Vc::frexp

  Double_v fexponent = IntToDouble(exponent-1-fNmin); //Backend
  //note: the normal  conversion from int to double,
  //Double_v fexponent(exponent-1-fNmin)
  //does not work for the int output of frexp which is in [int,dummy,int,dummy]

  ibin = Floor((mantissa-.5)*(2.*fNdiv)) + fNdiv*fexponent;

  Index_t idiv = ibin - fNdiv*Floor(ibin/fNdiv); //idiv = ibin%fNdiv
  //note: ibin%fNdiv = ibin & (fNdiv-1) for any fNdiv = 2^n does not work here
  //as the & operator is not vectorized)

  Double_v  power2  = ldexp(1.,exponent -1);
  Double_v  binsize = power2/fNdiv;
  Double_v  binloc  = power2 + binsize*idiv;

  frac = (x-binloc)/binsize;
}

} // end namespace impl
} // end namespace vecphys

#endif
