#ifndef GUComptonKleinNishina_H
#define GUComptonKleinNishina_H 1

#include "backend/Backend.h"

#include "GUConstants.h"
#include "GUTrack.h"

namespace vecphys {

VECPHYS_DEVICE_FORWARD_DECLARE( class GUAliasSampler; )
VECPHYS_DEVICE_DECLARE_CONV( GUAliasSampler )

inline namespace VECPHYS_IMPL_NAMESPACE {

class GUComptonKleinNishina
{
public:

  VECPHYS_CUDA_HEADER_BOTH
  GUComptonKleinNishina(Random_t* states, int threadId = -1); 

  VECPHYS_CUDA_HEADER_BOTH
  ~GUComptonKleinNishina();
  
  VECPHYS_CUDA_HEADER_BOTH
  void GetSampleParameters(double x, int &irow, int &icol,double &t);

  // Generate secondaries 

  VECPHYS_CUDA_HEADER_BOTH
  void Interact(GUTrack& projectile,    // In/Out: Updated to new state - choice
                int      targetElement, // Q: Need Material index instead ? 
                GUTrack* secondary )  const;

  VECPHYS_CUDA_HEADER_BOTH
  void InteractG4(GUTrack& inProjectile,
                  int      targetElement,
                  GUTrack* outSecondary );

#ifndef __CUDA_ARCH__  //  Method relevant only for *Vector* implementation
  // Vector version - stage 2 - Need the definition of these vector types
  void 
  Interact( // int ntracks,
            GUTrack_v& inProjectile,    // In/Out
            const int *targetElements,  // Number equal to num of tracks
            GUTrack_v* outSecondaryV    // Empty vector for secondaries
          ) const;     
          // TODO: use SOA arrays for output
#endif

  // Initializes this class and its sampler 
  VECPHYS_CUDA_HEADER_BOTH
  void BuildTable( int Z,
                              const double xmin,
                              const double xmax,
                              const int nrow,
			      const int ncol);
  // QUESTION: This might depend on physics? So maybe we should place inside model? 

  VECPHYS_CUDA_HEADER_BOTH
  void BuildPdfTable(int Z,
                     const double xmin,
                     const double xmax,
                     const int nrow,
                     const int ncol,
                     double *p);

  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH void 
  InteractKernel(typename Backend::Double_t energyIn, 
                 typename Backend::Double_t& energyOut,
		 typename Backend::Double_t& sinTheta) const;


  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 


  template<class Backend>
  VECPHYS_CUDA_HEADER_BOTH
  void SampleByCompositionRejection(typename Backend::Double_t energyIn,
				    typename Backend::Double_t& energyOut,
				    typename Backend::Double_t& sinTheta);

private: 
  // Implementation methods 

  VECPHYS_CUDA_HEADER_BOTH
  double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

  VECPHYS_CUDA_HEADER_BOTH
  void ConvertXtoFinalState(/* double energy, double sampledEphot, 
                               GUFourVector *newMomentum, GXTrack* secondary 
                             */) const;
  // Given the sampled output photon energy, create the output state




private:
  GUAliasSampler *fAliasSampler; 

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  Random_t* fRandomState;
  int       fThreadId;

  Precision fMinX;
  Precision fMaxX;
  Precision fDeltaX;

  Precision fMinY;
  Precision fMaxY;
  Precision fDeltaY;

  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void 
GUComptonKleinNishina::InteractKernel(typename Backend::Double_t energyIn, 
                                      typename Backend::Double_t& energyOut,
				      typename Backend::Double_t& sinTheta)
                                      const
{
  typedef typename Backend::Index_t  Index_t;
  typedef typename Backend::Double_t Double_t;

  Index_t   index;
  Index_t   icol;
  Double_t  fraction;

  fAliasSampler->SampleBin<Backend>(energyIn,index,icol,fraction);

  Double_t probNA ;
  Double_t aliasInd;

  //this does not work
  //  fAliasSampler->GatherAliasTable<Backend>(index,probNA,aliasInd);
  
  Double_t deltaE;
  deltaE = energyIn - energyIn/(1+2.0*energyIn*inv_electron_mass_c2);

  energyOut = fAliasSampler->SampleX<Backend>(deltaE,probNA,aliasInd,
  	 					       icol,fraction);
  sinTheta = SampleSinTheta<Backend>(energyIn,energyOut);
}    

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH
typename Backend::Double_t
GUComptonKleinNishina::
SampleSinTheta(typename Backend::Double_t energyIn,
               typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred photon

  Double_t epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > Backend::kOne;

  MaskedAssign( condition, 1.0 , &epsilon );

  Double_t E0_m    = inv_electron_mass_c2*energyIn;
  Double_t onecost = (Backend::kOne - epsilon)/(epsilon*E0_m);
  Double_t sint2   = onecost*(2.-onecost);

  Double_t sinTheta = 0.5;
  Bool_t condition2 = sint2 < Backend::kZero;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

template<class Backend>
VECPHYS_CUDA_HEADER_BOTH void
GUComptonKleinNishina::SampleByCompositionRejection(typename Backend::Double_t  energyIn,
						    typename Backend::Double_t& energyOut,
						    typename Backend::Double_t& sinTheta)
{
  typedef typename Backend::Double_t Double_t;

  Double_t epsilon, epsilonsq, onecost, sint2, greject ;
  
  Double_t E0_m = energyIn*inv_electron_mass_c2 ;
  Double_t eps0       = 1./(1. + 2.*E0_m);
  Double_t epsilon0sq = eps0*eps0;
  Double_t alpha1     = - log(eps0);
  Double_t alpha2     = 0.5*(1.- epsilon0sq);
  
  do {
    if( alpha1/(alpha1+alpha2) > UniformRandom(fRandomState,fThreadId) ) {
      epsilon   = exp(-alpha1*UniformRandom(fRandomState,fThreadId));
      epsilonsq = epsilon*epsilon; 
    } 
    else {
      epsilonsq = epsilon0sq+(1.- epsilon0sq)*UniformRandom(fRandomState,fThreadId);
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);
    
    //  } while (greject < 0.5);
  } while (greject < UniformRandom(fRandomState,fThreadId));
  
  energyOut = epsilon*energyIn;
  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

} // end namespace impl
} // end namespace vecphys

#endif
