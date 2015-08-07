#ifndef POWER2DIVISOR_H
#define POWER2DIVISOR_H 1

#include "backend/Backend.h"

namespace vecphys {

inline namespace VECPHYS_IMPL_NAMESPACE {

class Power2Divisor
{ 
public:
  VECPHYS_CUDA_HEADER_BOTH
  Power2Divisor(int nbins, int ndivs, Precision scale);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Index_t
  GetBin(typename Backend::Double_t x);

  template <typename Backend>
  inline 
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  GetBinPosition(typename Backend::Index_t ibin);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  GetBinSize(typename Backend::Index_t ibin);

  template <typename Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  FractionWithinBin(typename Backend::Double_t x,
                    typename Backend::Index_t ibin);

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetLowerBound() { return fLowerBound; }

  VECPHYS_CUDA_HEADER_BOTH
  Precision GetUpperBound() { return fUpperBound; }

private:
  int fNbins;              // number of bins
  int fNdivs;              // number of equal grids between 2^e and 2^{e+1}
  Precision fScale;        // scale for unit conversion (power of 10)

  Precision fLowerBound;   // lower bound of the range
  Precision fUpperBound;   // upper bound of the range

  Precision* fBinPosition; // bin position array
} ;

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Index_t
Power2Divisor::GetBin(typename Backend::Double_t x) 
{
  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Double_t Double_t;

  Int_t    exponent;
  Double_t mantissa = frexp (x/fScale, &exponent); // Vc::frexp
  Double_t fexponent = IntToDouble(exponent);
  return Floor((mantissa-.5)*(2.*fNdivs)) + 1. + fNdivs*(fexponent-1.);   
}

template <typename Backend>
inline
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
Power2Divisor::GetBinPosition(typename Backend::Index_t ibin) 
{
  //check for bounds
  return fBinPosition[ibin-1];   

  /* 
  // vectorized version - still expensive
  typedef typename Backend::Int_t Int_t;
  typedef typename Backend::Double_t Double_t;

  //no Vc % op, so do a%n = a - (n*int(a/n))
  Int_t exponent2 = DoubleToInt((ibin-1.)/fNdivs);
  Double_t mod = (ibin-1.)-fNdivs*IntToDouble(exponent2);

  //bin position = pow(2,int((i-1)/fNdivs))*(1.+((i-1)%fNdivs)/fNdivs) 
  Double_t binPosition = ldexp(1.+ mod/fNdivs, exponent2); // Vc::ldexp 
  return fScale*binPosition;
  */
}

// Specialisation for Vc Backend
#ifndef VECPHYS_NVCC
template<>
inline
VECPHYS_CUDA_HEADER_BOTH
kVc::Double_t
Power2Divisor::GetBinPosition<kVc>(typename kVc::Index_t ibin) 
{
  kVc::Double_t pos;
  for(int i = 0; i < kVc::kSize ; ++i) {
    int idx = ibin[i];
    pos[i] = fBinPosition[idx];
  }
  return pos;   
}
#endif

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
Power2Divisor::GetBinSize(typename Backend::Index_t ibin) 
{
  return GetBinPosition<Backend>(ibin+1) - GetBinPosition<Backend>(ibin);
}

template <typename Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
Power2Divisor::FractionWithinBin(typename Backend::Double_t x, 
                                 typename Backend::Index_t ibin) 
{
  typedef typename Backend::Double_t Double_t;

  Double_t down = GetBinPosition<Backend>(ibin);
  Double_t up   = GetBinPosition<Backend>(ibin+1);

  return (x - down)/(up-down);
}

} // end namespace impl
} // end namespace vecphys

#endif
