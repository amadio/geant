// 
//  Alias Sampler template function implementation
//   For use in implementing methods for scalar, vector, GPU
//   Depends on Backend 'technology' of VecGeom
//
//  First version using 'Backend' - 22 Oct 2014
//   Authors:  Sandro Wenzel, Soon Y. Jun, John Apostolakis
//
//  Desing choice: One Sampler per element 
//                    (potentially later per composit material?)
//
//  First version assumes linear 'x' integrand
//   TODO: a measure will be required for the integration for log, theta etc.
// 
//  Based on first alias sampler by Soon Y. Jun - July 2014

class GUAliasSampler
{
public: 
  GUAliasSampler(int    Zelement, 
                 double incomingMin, 
                 double incomingMax,
                 int    numEntriesIncoming,  // for 'energy' (or log) of projectile
                 int    numEntriesSampled 
                 );  

  template<class Backend>
  typename Backend::Double_t
  Sample( typename Backend::Double_t energyIn, 
          typename Backend::Double_t deltaY    ) const;  

  /* Builds (part) of our tables ( for one ene)
   *
  */
   // Builds all our table - must be called during initialisation 
  void SetPdf( double const * pdf );
  void BuildRow( int numVal, double* differentialXSection, double xMin, double xMax );

private:
// Implementation methods: 
  void BuildAliasTable( /*input: pdf constructed by Model */); // seems to be independent on model

  template<class Backend>
  typename Backend::Double_t
  Sample( typename Backend::Index_t  irow,   // ~ sampled value 
          typename Backend::Index_t  icol,   // ~ input Energy
          typename Backend::Double_t binSize,
          typename Backend::Double_t remainderX  //  in sampled variable
          ) const;

private:
  const int      fZelement; 
  
  const double   fIncomingMin; // Min of Incoming - e.g. e_Kinetic or Log(E_kinetic)
  const double   fIncomingMax; // Max
  const int      fInNumEntries;
  const double   fInverseBinIncoming;
  // const double   fInBinSize; 
  
  //  For the sampled variable 
  const int      fSampledNumEntries;   //  Old name fNcol  (number of Columns)
  const double   fSampledBinSize; 
  const double   fInverseBinSampled; 
  // Values needed per bin 
  //  -- Initialised how ??
  double*  fSampledMin; // Minimum value of 'x' sampled variable
  double*  fSampledMax; // Maximum value of 'x' sampled variable
  
  // Effective 2-dimensional arrays - size is fInNumEntries * fSampledNumEntries
  double * fpdf; // Original distribution
  
  double * fProbQ; // Non-alias probability ( sometimes called q in other code )
  int    * fAlias; // Alias table           -- could become Index_t ?
};

template<class Backend>
GUAliasSampler::
GUAliasSampler(int    Zelement, 
               double incomingMin, 
               double incomingMax,
               int    numEntriesIncoming,  // for 'energy' (or log) of projectile
               int    numEntriesSampled   
)  
  :
  fIncomingMin( incomingMin ),
  fIncomingMax( incomingMax ),
  fInNumEntries(numEntriesIncoming), 
  fInverseBinIncoming( numEntriesIncoming / incomingMax-incomingMin)
  fSampledNumEntries( numEntriesSampled ),
  fInverseBinSampled( 1.0 / (numEntriesSampled-1) )  // Careful - convention build / use table!
{
}

template<class Backend>
typename Backend::Double_t
GUAliasSampler::
Sample( typename Backend::Double_t energyIn ) const  
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Int_t    Int_t;	

  Int_t     irow, icol;
  Double_t  fraction;

  GetBin(energyIn, &irow, &icol, &fraction);

  // MinY = kineticEnergy/(1+2.0*kineticEnergy/electron_mass_c2);
  // MaxY = kineticEnergy;
  // DeltaY = fMaxY - fMinY;

  Double_t deltaY= fSampledMax[irow] - fSampledMin[irow]; 

  Double_t x = SampleX( irow, icol, fraction);

  return x; 
}

template<class Backend>
FQUALIFIER void
GUAliasSampler::GetBin(typename Backend::Double_t  kineticEnergy,
		       typename Backend::Int_t    &irow, 
		       typename Backend::Int_t    &icol,
		       typename Backend::Double_t &t) 
{
  irow = floor((kineticEnergy - fIncomingMin)*fInverseBinIncoming);
  G4double r1 = (fSampledNumEntries-1)*GPUniformRand(fDevState, fThreadId);
  icol = floor(r1);
  t = r1 - 1.0*icol;
}

template<class Backend>
typename Backend::Double_t
GUAliasSampler::SampleX( typename Backend::Index_t  irow, // Energy bin of incoming particle
                         typename Backend::Index_t  icol, // Bin index for this 'x' sample
                         typename Backend::Double_t rangeSampled, 
                         typename Backend::Double_t fraction )
{
  typedef typename Backend::Double_t Double_t;
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Int_t    Int_t;
  
  Double_t r1 = UniformRand(  );
  
  Double_t xd, xu;
  Double_t binSampled = rangeSampled * fInverseBinSampled; 
        // Was rangeSampled /(fSampledNumEntries-1);
  
  Double_t probNA;   // Non-alias probability
  Double_t aliasInd; //  This is really an integer -- could be Index_t !?  
  // fill probNA from table
  Index_t index = irow*fSampledNumEntries  + icol; 
  
  // Gather
  for( int i=0;i < Double_t::Size; ++i )
  {
    int iEntry= ConverfractionoInteger(index[i]);
    probNA[i]=    fPDFY[ iEntry ]; // index[i] ];
    aliasInd[i]=  fPDFA[ iEntry ]; // index[i] ];
  }
  // should investigate here whether gather is supported in Vc
  
  Bool_t condition = r1 <= probNA;
  
  // if branch
  
  MaskedAssign( condition, icol*binSampled ,     xd );   // Stores into xd
  MaskedAssign( condition, (icol+1)*binSampled , xu );   //        into xu
  
  // else branch
  
  MaskedAssign( !condition,  aliasInd*binSampled    , xd );
  MaskedAssign( !condition, (aliasInd+1)*binSampled , xu );
  
  Double_t x = (1 - fraction) * xd + fraction * xu;
  
  return x;
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
  double *ap = (double*)malloc(ncol*sizeof(double)); 

  //likelihood per equal probable event
  const double cp = 1.0/(ncol-1);

  for(int ir = 0; ir < nrow ; ++ir) {

    //copy and initialize
    for(int i = 0; i < nrow ; ++i) {
      fpdf[ir*ncol+i] = pdf[ir*ncol+i];
      a[i] = -1;
      ap[i] = p[i];
    }
  
    //O(n) iterations
    int iter = n;
  
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
      fpdfA[ir*ncol+recip] = donor;
      fpdfQ[ir*ncol+recip] = ncol*ap[recip];
    
      //update pdf 
      ap[donor] = ap[donor] - (cp-ap[recip]);
      ap[recip] = 0.0;

      --iter;
    }
    while (iter > 0);
  }

  free(ap);
}
