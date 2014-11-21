#include "GUAliasSampler.h"

GUAliasSampler::
GUAliasSampler(int    Zelement, 
               double incomingMin, 
               double incomingMax,
               int    numEntriesIncoming,  // for 'energy' (or log) of projectile
               int    numEntriesSampled   
)  
  :
  fZelement(Zelement),
  fIncomingMin( incomingMin ),
  fIncomingMax( incomingMax ),
  fInNumEntries(numEntriesIncoming), 
  fInverseBinIncoming( numEntriesIncoming / (incomingMax-incomingMin)),
  fSampledNumEntries( numEntriesSampled ),
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) ),  // Careful - convention build / use table!
  fSampledBinSize(1.0/* WHAT TO PUT HERE ?? */)
{
  fpdf = new double [fInNumEntries*fSampledNumEntries];
  fProbQ = new double [fInNumEntries*fSampledNumEntries];
  fAlias = new int [fInNumEntries*fSampledNumEntries];
}

GUAliasSampler::~GUAliasSampler()
{
  if(fpdf)   delete [] fpdf;
  if(fProbQ) delete [] fProbQ;
  if(fAlias) delete [] fAlias;
}

void GUAliasSampler::BuildAliasTables( const int nrow,
                                       const int ncol,
                                       double   *pdf )
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
  int *a = (int*)malloc(ncol*sizeof(int)); 
  double *ap = (double*)malloc(ncol*sizeof(double)); 

  //likelihood per equal probable event
  const double cp = 1.0/(ncol-1);

  for(int ir = 0; ir < nrow ; ++ir) {

    //copy and initialize
    for(int i = 0; i < ncol ; ++i) {
      fpdf[ir*ncol+i] = pdf[ir*ncol+i];
      a[i] = -1;
      ap[i] = pdf[ir*ncol+i];
    }

    //O(n) iterations
    int iter = ncol;
  
    do {
      int donor = 0;
      int recip = 0;
    
      //have the better search algorithm?
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
      fAlias[ir*ncol+recip] = donor;
      fProbQ[ir*ncol+recip] = ncol*ap[recip];
    
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
