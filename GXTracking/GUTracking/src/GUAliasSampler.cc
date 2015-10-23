#include "GUAliasSampler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
GUAliasSampler::
GUAliasSampler(Random_t* states, int threadId,
               int    maxZelement, 
               double incomingMin, 
               double incomingMax,
               int    numEntriesIncoming, // for 'energy' (or log) of projectile
               int    numEntriesSampled   
)  
  :
  fRandomState(states), 
  fThreadId(threadId),
  fMaxZelement(maxZelement),
  fIncomingMin( incomingMin ),
  fIncomingMax( incomingMax ),
  fInNumEntries(numEntriesIncoming), 
  fInverseBinIncoming( numEntriesIncoming / (incomingMax-incomingMin)),
  fInverseLogBinIncoming( numEntriesIncoming / (log(incomingMax)-log(incomingMin))),
  fSampledNumEntries( numEntriesSampled ),
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) )
{

  fAliasTableManager = new GUAliasTableManager(maxZelement);

}

VECPHYS_CUDA_HEADER_BOTH
GUAliasSampler::
GUAliasSampler(Random_t* states, int threadId,
               // int    Zelement, 
               double incomingMin, 
               double incomingMax,
               int    numEntriesIncoming, // for 'energy' (or log) of projectile
               int    numEntriesSampled,   
               GUAliasTableManager* tableManager
)  
  :
  fRandomState(states), 
  fThreadId(threadId),
  // fZelement(Zelement),
  fIncomingMin( incomingMin ),
  fIncomingMax( incomingMax ),
  fInNumEntries(numEntriesIncoming), 
  fInverseBinIncoming( numEntriesIncoming / (incomingMax-incomingMin)),
  fInverseLogBinIncoming( numEntriesIncoming / (log(incomingMax)-log(incomingMin))),
  fSampledNumEntries( numEntriesSampled ),
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) )  // Careful - convention build / use table!
{
  fAliasTableManager = tableManager;
}

VECPHYS_CUDA_HEADER_BOTH
GUAliasSampler::~GUAliasSampler()
{
#ifndef VECPHYS_NVCC
  if(fAliasTableManager)  delete fAliasTableManager;
#endif

}

VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::PrintTable()
{
  printf("Incoming Min= %g , Max= %g , numEntries= %d \n",
         fIncomingMin, fIncomingMax, fInNumEntries );

  if( fAliasTableManager->GetNumberOfElements() >0 ) {
    for (int i = 0; i < fAliasTableManager->GetNumberOfElements() ; ++i ) {
      GUAliasTable* tb= fAliasTableManager->GetAliasTable(i);
      printf("GUAliasSampler fAliasTable = %p fNGrid = %d value=(%d,%f,%f)\n",
	     tb,
	     tb->fNGrid,
	     tb->fAlias[1],
	     tb->fProbQ[1],
	     tb->fpdf[1]);
    }
  }
  else {
    printf("GUAliasSampler fAliasTableManager is empty\n");
  }

}
VECPHYS_CUDA_HEADER_HOST
void GUAliasSampler::BuildAliasTable( int Zelement,
                                      int nrow,
                                      int ncol,
                                      const double *pdf )
{
  // Build alias and alias probability
  //    
  // Reference: (1) A.J. Walker, "An Efficient Method for Generating Discrete 
  // Random Variables with General Distributions" ACM Trans. Math. Software, 3,
  // 3, 253-256 (1977) (2) A.L. Edwards, J.A. Rathkopf, and R.K. Smidt, 
  // "Extending the Alias Monte Carlo Sampling Method to General Distributions"
  // UCRL-JC-104791 (1991)
  //
  // input :  nrow             (multiplicity of alias tables)
  //          ncol             (dimension of discrete outcomes)
  //          pdf[nrow x ncol] (probability density function)
  // output:  a[nrow x ncol]   (alias)
  //          q[nrow x ncol]   (non-alias probability) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //temporary array
  int *a     = (int*)   malloc(ncol*sizeof(int)); 
  double *ap = (double*)malloc(ncol*sizeof(double)); 

  //likelihood per equal probable event
  const double cp = 1.0/(ncol-1);

  GUAliasTable* table = new GUAliasTable((nrow+1)*ncol);

  for(int ir = 0; ir <= nrow ; ++ir) {

    //copy and initialize
    for(int i = 0; i < ncol ; ++i) {

       int ipos= ir*ncol+i;
       table->fpdf[ipos] = pdf[ipos];

       a[i] = -1;
       ap[i] = pdf[ipos]; // pdf[ir*ncol+i];
    }

    //O(n) iterations
    int iter = ncol;
  
    do {
      int donor = 0;
      int recip = 0;
    
      // A very simple search algorithm
      for(int j = donor; j < ncol ; ++j) {
         if(ap[j] >= cp) {
            donor = j;
            break;
         }
      }

      for(int j = recip; j < ncol ; ++j) {
         if(ap[j] > 0.0 && ap[j] < cp) {
            recip = j;
            break;
         }
      }

      //alias and non-alias probability

      table->fAlias[ir*ncol+recip] = donor;
      table->fProbQ[ir*ncol+recip] = ncol*ap[recip];
    
      //update pdf 
      ap[donor] = ap[donor] - (cp-ap[recip]);
      ap[recip] = 0.0;
      --iter;

    }
    while (iter > 0);
  }

  fAliasTableManager->AddAliasTable(Zelement,table);

  free(a);
  free(ap);
}

} // end namespace impl
} // end namespace vecphys
