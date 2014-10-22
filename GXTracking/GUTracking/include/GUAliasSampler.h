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
//  Based on first alias sampler by Soon Y. Jun - July 2014

class GUAliasSampler
{
  public: 
  	GUAliasSampler(int Zelement);  // One table per element !!

    template<class Backend>
    typename Backend::Double_t
    SampleX( typename Backend::Int_t irow,   // ~ sampled value 
		     typename Backend::Int_t icol,   // ~ input Energy
			 typename Backend::Double_t dy, 
			 typename Backend::Double_t tt );
    BuildTable();  
      // Builds all our table - must be called during initialisation 

private:
// Implementation methods: 
	BuildPdfTable(); 
	BuildAlisTable(); 

private:
    int      fZelement; 
	double * fpdfX; // Original distribution
	double * fpdfY; // Non-alias probability 
	int    * fpdfA; // Alias table
};


template<class Backend>
typename Backend::Double_t
GUAliasSampler::SampleX( typename Backend::Int_t irow, 
	   typename Backend::Int_t icol, 
	   typename Backend::Double_t dy, 
	   typename Backend::Double_t tt )
{
	typedef typename Backend::Double_t Double_t;
	typedef typename Backend::Bool_t Bool_t;
	typedef typename Backend::Int_t Int_t;
	
	
	Double_t r1 = UniformRand(  );
	
	Double_t xd, xu;
	Double_t dx = dy/(fNcol-1); // dy*(fNcolMinusInverse)

	Double_t thresh;
	Double_t Aval;
	// fill thresh from table
	Int_t index = irow*fNcol  + icol; 
	
	// should investigate here whether gather is supported in Vc
	for( int i=0;i < Double_t::Size; ++i )
	{
		thresh[i]=fPDFY[ index[i] ];
	    Aval[i]=fPDFA[ index[i] ];
	}
	
	Bool_t condition = r1 <= thresh;
	
	// if branch
	
	MaskedAssign( condition, dx*icol ,xd );
	MaskedAssign( condition, dx*(icol+1) ,xu );

	// else branch
	
	MaskedAssign( !condition, dx*Aval ,xd );
    MaskedAssign( !condition, dx*(Aval + 1) ,xu );

    Double_t x = ( 1 - tt)*xd + tt*xu;

    return x;
}
