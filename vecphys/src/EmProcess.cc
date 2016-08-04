#include "EmProcess.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

// scalar implementation without Backend

template <class Process>
VECCORE_CUDA_HOST void EmProcess<Process>::BuildAlias()
{
  // build alias table for the physics process based on their relative cross sections
  // non-alias proability will be calculated on the fly

  for (int i = 0; i < fNumberOfMaterialBin; ++i) {
    for (int j = 0; j < fNumberOfEnergyBin; ++j) {

      int ibin = i * fNumberOfEnergyBin + j;

      // alias table
      int *a = new int[fNumberOfProcess];
      double *ap = new double[fNumberOfProcess];

      const double cp = 1. / fNumberOfProcess;

      // copy and initialize
      double pdf[3] = {fCrossSectionData[ibin].fWeight[0], fCrossSectionData[ibin].fWeight[1],
                       1.0 - fCrossSectionData[ibin].fWeight[0] - fCrossSectionData[ibin].fWeight[1]};

      for (int k = 0; k < fNumberOfProcess; ++k) {
        a[k] = -1;
        ap[k] = pdf[k];
      }
      // O(n) iterations
      int iter = fNumberOfProcess;

      do {
        int donor = 0;
        int recip = 0;

        // A very simple search algorithm
        for (int k = donor; k < fNumberOfProcess; ++k) {
          if (ap[k] >= cp) {
            donor = k;
            break;
          }
        }

        for (int k = recip; k < fNumberOfProcess; ++k) {
          if (ap[k] >= 0.0 && ap[k] < cp) {
            recip = k;
            break;
          }
        }

        // alias and non-alias probability
        fCrossSectionData[ibin].fAlias[recip] = donor;

        // update pdf
        ap[donor] = ap[donor] - (cp - ap[recip]);
        ap[recip] = -1.0;
        --iter;

      } while (iter > 0);

      delete [] a;
      delete [] ap;
    }
  }
}

} // end namespace impl
} // end namespace vecphys
