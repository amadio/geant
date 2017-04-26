#include "GUAliasSampler.h"

#include "MaterialHandler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
GUAliasSampler::GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                               int numEntriesIncoming, // for 'energy' (or log) of projectile
                               int numEntriesSampled)
    : fRandomState(states), fThreadId(threadId), fIncomingMin(incomingMin), fIncomingMax(incomingMax),
      fInNumEntries(numEntriesIncoming), fLogIncomingMin(math::Log(incomingMin)),
      fInverseBinIncoming(numEntriesIncoming / (incomingMax - incomingMin)),
      fInverseLogBinIncoming(numEntriesIncoming / (math::Log(incomingMax) - fLogIncomingMin)),
      fSampledNumEntries(numEntriesSampled), fInverseBinSampled(1.0 / numEntriesSampled)
{
  int nelements = MaterialHandler::Instance()->GetNumberOfElements();
  int ngrid = (fInNumEntries + 1) * fSampledNumEntries;
  fAliasTableManager = new GUAliasTableManager(nelements, ngrid);
}

VECCORE_ATT_HOST_DEVICE
GUAliasSampler::GUAliasSampler(Random_t *states, int threadId, double incomingMin, double incomingMax,
                               int numEntriesIncoming, // for 'energy' (or log) of projectile
                               int numEntriesSampled, GUAliasTableManager *tableManager)
    : fRandomState(states), fThreadId(threadId), fIncomingMin(incomingMin), fIncomingMax(incomingMax),
      fInNumEntries(numEntriesIncoming), fLogIncomingMin(math::Log(incomingMin)),
      fInverseBinIncoming(numEntriesIncoming / (incomingMax - incomingMin)),
      fInverseLogBinIncoming(numEntriesIncoming / (math::Log(incomingMax) - fLogIncomingMin)),
      fSampledNumEntries(numEntriesSampled),
      fInverseBinSampled(1.0 / numEntriesSampled) // Careful - convention build / use table!
{
  fAliasTableManager = tableManager;
}

VECCORE_ATT_HOST_DEVICE
GUAliasSampler::~GUAliasSampler()
{
#if !defined(VECCORE_CUDA)
  if (fAliasTableManager)
    delete fAliasTableManager;
#endif
}

VECCORE_ATT_HOST_DEVICE
void GUAliasSampler::PrintTable()
{
  printf("Incoming Min= %g , Max= %g , numEntries= %d \n", fIncomingMin, fIncomingMax, fInNumEntries);

  if (fAliasTableManager->GetNumberOfElements() > 0) {
    for (int i = 0; i < fAliasTableManager->GetNumberOfElements(); ++i) {
      GUAliasTable *tb = fAliasTableManager->GetAliasTable(i);
      printf("GUAliasSampler fAliasTable = %p fNGrid = %d value=(%d,%f,%f)\n", tb, tb->fNGrid, tb->fAlias[1],
             tb->fProbQ[1], tb->fpdf[1]);
    }
  }
  else {
    printf("GUAliasSampler fAliasTableManager is empty\n");
  }
}
VECCORE_ATT_HOST
void GUAliasSampler::BuildAliasTable(int Zelement, const double *pdf)
{
  // Build alias and alias probability
  //
  // Reference: (1) A.J. Walker, "An Efficient Method for Generating Discrete
  // Random Variables with General Distributions" ACM Trans. Math. Software, 3,
  // 3, 253-256 (1977) (2) A.L. Edwards, J.A. Rathkopf, and R.K. Smidt,
  // "Extending the Alias Monte Carlo Sampling Method to General Distributions"
  // UCRL-JC-104791 (1991)
  //
  // input : fInNumEntries       multiplicity of alias tables)
  //         fSampledNumEntries (dimension of discrete outcomes)
  //         pdf[fInNumEntries x fSampledNumEntries] (probability density function)
  // output: a[fInNumEntries x fSampledNumEntries]   (alias)
  //         q[fInNumEntries x fSampledNumEntries]   (non-alias probability)
  //

  // temporary array
  double *ap = (double *)malloc(fSampledNumEntries * sizeof(double));

  // likelihood per equal probable event
  const double cp = 1.0 / fSampledNumEntries;

  GUAliasTable *table = new GUAliasTable((fInNumEntries + 1) * fSampledNumEntries);

  for (int ir = 0; ir <= fInNumEntries; ++ir) {

    // copy and initialize
    for (int i = 0; i < fSampledNumEntries ; ++i) {

      int ipos = ir * fSampledNumEntries + i;
      table->fpdf[ipos] = pdf[ipos];

      ap[i] = pdf[ipos]; // pdf[ir*fSampledNumEntries+i];
    }

    // O(n) iterations
    int iter = fSampledNumEntries;

    do {
      int donor = 0;
      int recip = 0;

      // A very simple search algorithm
      for (int j = donor; j < fSampledNumEntries; ++j) {
        if (ap[j] >= cp) {
          donor = j;
          break;
        }
      }

      for (int j = recip; j < fSampledNumEntries; ++j) {
        if (ap[j] > 0.0 && ap[j] < cp) {
          recip = j;
          break;
        }
      }

      // alias and non-alias probability

      table->fAlias[ir * fSampledNumEntries + recip] = donor;
      table->fProbQ[ir * fSampledNumEntries + recip] = fSampledNumEntries * ap[recip];

      // update pdf
      ap[donor] = ap[donor] - (cp - ap[recip]);
      ap[recip] = 0.0;
      --iter;

    } while (iter > 0);
  }

  fAliasTableManager->AddAliasTable(Zelement, table);

  free(ap);
}

} // end namespace impl
} // end namespace vecphys
