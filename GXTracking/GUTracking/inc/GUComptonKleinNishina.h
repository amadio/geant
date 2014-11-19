#ifndef GUComptonKleinNishina_H
#define GUComptonKleinNishina_H 1

#include "GUTypeDef.h"
#include "GUTrack.h"

class GUAliasSampler; 

class GUComptonKleinNishina
{
public:

  FQUALIFIER GUComptonKleinNishina(); 
  FQUALIFIER ~GUComptonKleinNishina();
  
  FQUALIFIER void GetSampleParameters(double x, int &irow, int &icol,
				      double &t);
  // Generate secondaries 
  FQUALIFIER 
  void Interact(GUTrack& projectile,    // In/Out: Updated to new state - choice
                int      targetElement, // Q: Need Material index instead ? Both int => careful
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
  // configure the 'mode' ??
  // FQUALIFIER void UseOriginalSampler(bool useOriginal);  // 
  // FQUALIFIER bool UsingOriginalSampler() const; 

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

  FQUALIFIER double SampleDirection(double energyIn, double energyOutPhoton);

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

#endif

