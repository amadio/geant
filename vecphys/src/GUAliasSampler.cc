#include "GUAliasSampler.h"

#include "MaterialHandler.h"
#include "Geant/Error.h"

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
    int iter = fSampledNumEntries + 1;

    do {
      int donor = 0;
      int recip = 0;

      // A very simple algorithm building the alias table

      // donor : note that ap[j] = cp is the case for a flat distribution
      for (int j = donor ; j < fSampledNumEntries; ++j) {
	if (ap[j] >= cp ) {
          donor = j;
          break;
        }
      }

      // recipeant
      for (int j = recip ; j < fSampledNumEntries; ++j) {
        if (ap[j] > 0. && ap[j] <= cp) {
          recip = j;
          break;
        }
      }

      // alias and non-alias probability
      table->fAlias[ir * fSampledNumEntries + recip] = donor;
      table->fProbQ[ir * fSampledNumEntries + recip] = fSampledNumEntries * ap[recip];

      // update pdf of the donor
      ap[donor] = ap[donor] - (cp - ap[recip]);

      // remove the recipeant from the search list      
      ap[recip] = 0.0;

      --iter;

    } while (iter > 0);

    // check the validity of the table
    for (int i = 0 ; i < fSampledNumEntries ; ++i) {
      if(table->fAlias[ir * fSampledNumEntries + i] < 0 || table->fProbQ[ir * fSampledNumEntries + i] < -1.0e-10) {
	// this algorithm has a flaw
	Geant::Error("GUAliasSampler::BuildAliasTable", "Invalid alias table entries");
	//        printf("GUAliasSampler::BuildAliasTable : ERROR building the alias table\n");
	//        printf("  (fInNumEntries,fSampledNumEntries)=(%d,%d)\n",ir,i);
	//        printf("  (fAlias,fProbQ)=(%d,%f)\n",table->fAlias[ir * fSampledNumEntries + i],
	//                                             table->fProbQ[ir * fSampledNumEntries + i]);
      }
    }
  }

  fAliasTableManager->AddAliasTable(Zelement, table);

  free(ap);
}

} // end namespace impl
} // end namespace vecphys
