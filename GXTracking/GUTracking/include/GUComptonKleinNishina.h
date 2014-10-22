#ifndef GUComptonKleinNishina_H
#define GUComptonKleinNishina_H 1

#include "GUTypeDef.h"
#include "GPRandom.h"
// #include "GPThreeVector.h"
#include "GXTrack.h"

// #include "GXVDiscreteSampling.h"
class GUAliasSampler; 

class GUComptonKleinNishina  // : public GXVDiscreteSampling
{
public:

  FQUALIFIER GUComptonKleinNishina(); 

  /*----------------
			    double minE,   // Minimum E of model & table
			    double maxE,   // Maximum E of model & table

			    double dimRow,    // Number of table entries for one E
			    int dimCol,    // Number of Energies in grid
			    double *pdfX,  // Table of Original PDF 
			    double *pdfY   // Table of Inverse Cumulative distribution

          curandState* devState,
          int threadId,

			    int    *pdfA,
			    double *pdfX,
			    double *pdfY);
    ---------------------
    */
  FQUALIFIER ~GXKleinNishina();
  
  FQUALIFIER void GetSampleParameters(double x, int &irow, int &icol,
				      double &t);
    // Generate secondaries 
  FQUALIFIER 
  void Interact(double   energy,  GUFourVector   *newMomentum,  GXTrack*   secondary)  const;

#ifndef __CUDA_ARCH__  //  Method relevant only for *Vector* implementation
  // void Interact(const double_v energyV, GUFourVector_v *newMomentumV, GXTrack_v* secondaryV) const;     
  // Vector version - stage 2 - Need the definition of these vector types
  void Interact(int ntracks, 
                const double* energyV, 
                std::vector<GUFourVector> &outMomentumV,            
                std::vector<GXTrack_v>    &outSecondaryV  ) const;     
               // TODO: use SOA arrays for output
#endif
  // configure the 'mode' ??
  // FQUALIFIER void UseOriginalSampler(bool useOriginal);  // 
  // FQUALIFIER bool UsingOriginalSampler() const; 

    // Initializes this class and its sampler 
  FQUALIFIER void BuildTables(); 
private: 
  // Implementation methods 
  FQUALIFIER
  void ConvertXtoFinalState(double energy, double sampledEphot, GUFourVector *newMomentum,  GXTrack*   secondary)  const;
   // Given the sampled output photon energy, create the output state

  FQUALIFIER double SampleDirection(double energyIn, double energyOutPhoton);

private:
  GUAliasSampler *fAliasSampler; 
  XSectionKleinNishina *fXSection;

  // Helper data members for GPU random -- to be replaced by use of a GPU manager class
  int fThreadId;
  // curandState* fDevState;

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

