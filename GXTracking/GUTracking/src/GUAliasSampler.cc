#include "GUAliasSampler.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_BOTH
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
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) )  // Careful - convention build / use table!
  // fSampledBinSize(1.0/* WHAT TO PUT HERE ?? */)
{
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  // For multiple elements, this should be a vector of GUAliasTable
  // fAliasTable = new (GUAliasTable*)[maxZelement+1];
  fAliasTable = (GUAliasTable**) malloc( (maxZelement+1)*sizeof(GUAliasTable*) ); 

  for( int z=1; z<maxZelement; z++)
  {
     // std::cout << " Creating fAliasTable for Z= " << z << std::endl;
     fAliasTable[z] = new GUAliasTable(fInNumEntries*fSampledNumEntries);
     // std::cout << " Pointer for Z= " << z << " = > " << fAliasTable[z] << std::endl;
  }
}

VECPHYS_CUDA_HEADER_BOTH
GUAliasSampler::
GUAliasSampler(Random_t* states, int threadId,
               // int    Zelement, 
               double incomingMin, 
               double incomingMax,
               int    numEntriesIncoming, // for 'energy' (or log) of projectile
               int    numEntriesSampled,   
               GUAliasTable* table
)  
  :
  fRandomState(states), 
  fThreadId(threadId),
  // fZelement(Zelement),
  fIncomingMin( incomingMin ),
  fIncomingMax( incomingMax ),
  fInNumEntries(numEntriesIncoming), 
  fInverseBinIncoming( numEntriesIncoming / (incomingMax-incomingMin)),
  fSampledNumEntries( numEntriesSampled ),
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) )  // Careful - convention build / use table!
  // fSampledBinSize(1.0 ) // ?? WHAT TO PUT HERE ?? 
{
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  // For multiple elements, this should be a vector of GUAliasTable
  fAliasTable[0] = table;
}

VECPHYS_CUDA_HEADER_BOTH
GUAliasSampler::~GUAliasSampler()
{
   if(fAliasTable)
   {
      for( int z=1; z< fMaxZelement ; ++z )
         delete fAliasTable[z];

      delete fAliasTable;
      fAliasTable= (GUAliasTable**) 0;
   }
}

VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::PrintTable()
{
  if(fAliasTable != NULL)
  {
     printf("GUAliasSampler fAliasTable = %p \n", fAliasTable ); 
     for ( int z=1; z< fMaxZelement; ++z )
     {
        GUAliasTable* tb= fAliasTable[z];
        printf("GUAliasSampler fAliasTable = %p fNGrid = %d value=(%d,%f,%f)\n",
               tb,
               tb->fNGrid,
               tb->fAlias[1],
               tb->fProbQ[1],
               tb->fpdf[1]);
     }
  }
  else {
    printf("GUAliasSampler fAliasTable is NULL\n");
  }
}

VECPHYS_CUDA_HEADER_BOTH
void GUAliasSampler::GetAlias(int    index,
                              int    zElement,
                              double &probNA,  
                              int    &aliasInd) const 
{
  assert( zElement > 0  && zElement <= fMaxZelement );
  
  //gather for alias table lookups
  probNA =    fAliasTable[zElement]->fProbQ[ index ];
  aliasInd =  fAliasTable[zElement]->fAlias[ index ];
}

VECPHYS_CUDA_HEADER_BOTH
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

  // std::cout << " GUAliasSampler - building table for element Z= " << Zelement << std::endl;
   
  //temporary array
  int *a     = (int*)   malloc(ncol*sizeof(int)); 
  double *ap = (double*)malloc(ncol*sizeof(double)); 

  //likelihood per equal probable event
  const double cp = 1.0/(ncol-1);

  // std::cout << " Pointer fAliasTable[ Z= " << Zelement << " ] > " << fAliasTable[Zelement] << std::endl; 
  // std::cout << "  fAliasTable[ Z= " << Zelement << " ]->fpdf = " << fAliasTable[Zelement]->fpdf << std::endl;

  GUAliasTable* aliasTable= fAliasTable[Zelement];

  for(int ir = 0; ir < nrow ; ++ir) {

    //copy and initialize
    for(int i = 0; i < ncol ; ++i) {

       // fAliasTable[Zelement]->fpdf[ir*ncol+i] = pdf[ir*ncol+i];
       int ipos= ir*ncol+i;
       aliasTable->fpdf[ipos] = pdf[ipos];
       // aliasTable->fpdf[ir*ncol+i] = pdf[ir*ncol+i];       

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

      // fAliasTable[Zelement]         
      aliasTable
         ->fAlias[ir*ncol+recip] = donor;
      // fAliasTable[Zelement]
      aliasTable
         ->fProbQ[ir*ncol+recip] = ncol*ap[recip];
    
      //update pdf 
      ap[donor] = ap[donor] - (cp-ap[recip]);
      ap[recip] = 0.0;
      --iter;

    }
    while (iter > 0);
  }

  free(a);
  free(ap);
}

} // end namespace impl
} // end namespace vecphys
