#include "Power2Divisor.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_BOTH
Power2Divisor::Power2Divisor(int nbins, int ndivs, Precision scale) 
  : fNbins(nbins), fNdivs(ndivs), fScale(scale)
{
  fBinPosition = new Precision [fNbins];

  for(int i = 0 ; i < fNbins ; ++i) {
    fBinPosition[i] = ldexp((1.+((i-1)%fNdivs)/(fNdivs)),int((i-1)/fNdivs));
  }

  fLowerBound = fScale;
  fUpperBound = fScale*fBinPosition[fNbins-1];//GetBinPosition<kScalar>(fNbins);
}

} // end namespace impl
} // end namespace vecphys
