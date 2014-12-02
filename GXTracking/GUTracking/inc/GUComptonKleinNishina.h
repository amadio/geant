#ifndef GUComptonKleinNishina_H
#define GUComptonKleinNishina_H 1

#include "GUTypeDef.h"
#include "GUConstants.h"
#include "GUTrack.h"
#include "GURandom.h"

class GUAliasSampler; 

class GUComptonKleinNishina
{
public:

  FQUALIFIER GUComptonKleinNishina(); 
  FQUALIFIER ~GUComptonKleinNishina();
  
  FQUALIFIER 
  void GetSampleParameters(double x, int &irow, int &icol,double &t);

  // Generate secondaries 
  FQUALIFIER 
  void Interact(GUTrack& projectile,    // In/Out: Updated to new state - choice
                int      targetElement, // Q: Need Material index instead ? 
                GUTrack* secondary )  const;

#ifndef __CUDA_ARCH__  //  Method relevant only for *Vector* implementation
  // Vector version - stage 2 - Need the definition of these vector types
  FQUALIFIER void 
  Interact( // int ntracks,
            GUTrack_v& inProjectile,    // In/Out
            const int *targetElements,  // Number equal to num of tracks
            GUTrack_v* outSecondaryV    // Empty vector for secondaries
          ) const;     
          // TODO: use SOA arrays for output
#endif
  FQUALIFIER void SampleByCompositionRejection(double energyIn,
					       double& energyOut,
					       double& sinTheta);
  // Initializes this class and its sampler 
  FQUALIFIER void BuildTable( int Z,
                              const double xmin,
                              const double xmax,
                              const int nrow,
			      const int ncol);
  // QUESTION: This might depend on physics? So maybe we should place inside model? 

  FQUALIFIER void BuildPdfTable(int Z,
                                const double xmin,
                                const double xmax,
                                const int nrow,
                                const int ncol,
                                double *p);

private: 
  // Implementation methods 
  FQUALIFIER double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

 FQUALIFIER
   void ConvertXtoFinalState(/* double energy, double sampledEphot, 
				GUFourVector *newMomentum, GXTrack* secondary 
                              */) const;
  // Given the sampled output photon energy, create the output state

  template<class Backend>
  typename Backend::Double_t
  SampleSinTheta(typename Backend::Double_t energyIn,
                 typename Backend::Double_t energyOut) const; 

private:
  GUAliasSampler *fAliasSampler; 
  //  GUXSectionKleinNishina *fXSection;

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  randomState* fRandomState;
  int fThreadId;

  double fMinX;
  double fMaxX;
  double fDeltaX;

  double fMinY;
  double fMaxY;
  double fDeltaY;

  //Sampling Tables
  int fNrow;
  int fNcol;
};

//Implementation

template<class Backend>
typename Backend::Double_t
GUComptonKleinNishina::
SampleSinTheta(typename Backend::Double_t energyIn,
               typename Backend::Double_t energyOut) const
{
  typedef typename Backend::Bool_t   Bool_t;
  typedef typename Backend::Double_t Double_t;

  //angle of the scatterred photon

  Double_t epsilon = energyOut/energyIn;

  Bool_t condition = epsilon > Vc::One;
  MaskedAssign( condition, 1.0 , &epsilon );

  Double_t E0_m    = inv_electron_mass_c2*energyIn;
  Double_t onecost = (Vc::One - epsilon)/(epsilon*E0_m);
  Double_t sint2   = onecost*(2.-onecost);

  Double_t sinTheta;
  Bool_t condition2 = sint2 < Vc::Zero;

  MaskedAssign(  condition2, 0.0, &sinTheta );   // Set sinTheta = 0
  MaskedAssign( !condition2, Vc::sqrt(sint2), &sinTheta );   
  
  return sinTheta;
}

#endif
