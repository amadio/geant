#ifndef POWER2DIVISOR_H
#define POWER2DIVISOR_H 1

#include "base/Global.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

class Power2Divisor
{
public:
  VECPHYS_CUDA_HEADER_BOTH
  Power2Divisor(int nmin, int nmax, int ndiv);

  VECPHYS_CUDA_HEADER_BOTH
  int GetNumberOfBins();

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetLowerBound() { return fLowerBound; }

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetUpperBound() { return fUpperBound; }

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetBinPosition(int ibin);

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetBinSize(int ibin);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void GetBinAndFraction(typename Backend::Double_t x,
                         typename Backend::Index_t& ibin,
                         typename Backend::Double_t& frac);

private:
  VECPHYS_CUDA_HEADER_BOTH
  int Power2Exponent(int ibin);

private:
  int fNmin;              // minimum of the power2 exponent
  int fNmax;              // maximum of the power2 exponent
  int fNdiv;              // number of equal divisions between 2^e and 2^{e+1}

  Precision fLowerBound;   // lower bound of the power2 range
  Precision fUpperBound;   // upper bound of the power2 range
};

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH
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
