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
  Precision GetLowerBound() { return fLowerBound; }

  VECCORE_CUDA_HOST_DEVICE
  Precision GetUpperBound() { return fUpperBound; }

  VECCORE_CUDA_HOST_DEVICE
  Precision GetBinPosition(int ibin);

  VECCORE_CUDA_HOST_DEVICE
  Precision GetBinSize(int ibin);

  template <typename Backend>
  VECCORE_CUDA_HOST_DEVICE
  void GetBinAndFraction(typename Backend::Double_t x,
                         typename Backend::Index_t& ibin,
                         typename Backend::Double_t& frac);

private:
  VECCORE_CUDA_HOST_DEVICE
  int Power2Exponent(int ibin);

private:
  int fNmin;              // minimum of the power2 exponent
  int fNmax;              // maximum of the power2 exponent
  int fNdiv;              // number of equal divisions between 2^e and 2^{e+1}

  Precision fLowerBound;   // lower bound of the power2 range
  Precision fUpperBound;   // upper bound of the power2 range
};

template <typename Backend>
VECCORE_CUDA_HOST_DEVICE
void Power2Divisor::GetBinAndFraction(typename Backend::Double_t x,
                                      typename Backend::Index_t& ibin,
                                      typename Backend::Double_t& frac)
{
  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Index_t Index_t;
  typedef typename Backend::Double_t Double_t;

  Int_t    exponent;
  Double_t mantissa = frexp (x, &exponent); // Vc::frexp

  Double_t fexponent = IntToDouble(exponent-1-fNmin); //Backend
  //note: the normal  conversion from int to double,
  //Double_t fexponent(exponent-1-fNmin)
  //does not work for the int output of frexp which is in [int,dummy,int,dummy]

  ibin = Floor((mantissa-.5)*(2.*fNdiv)) + fNdiv*fexponent;

  Index_t idiv = ibin - fNdiv*Floor(ibin/fNdiv); //idiv = ibin%fNdiv
  //note: ibin%fNdiv = ibin & (fNdiv-1) for any fNdiv = 2^n does not work here
  //as the & operator is not vectorized)

  Double_t  power2  = ldexp(1.,exponent -1);
  Double_t  binsize = power2/fNdiv;
  Double_t  binloc  = power2 + binsize*idiv;

  frac = (x-binloc)/binsize;
}

} // end namespace impl
} // end namespace vecphys

#endif
