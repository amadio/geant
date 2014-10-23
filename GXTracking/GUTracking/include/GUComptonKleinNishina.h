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
  FQUALIFIER void BuildPdfTable(); // QUESTION: This might depend on physics? So maybe we should place inside model? 
  FQUALIFIER void BuildTables(); 
private: 
  // Implementation methods 
  FQUALIFIER double CalculateDiffCrossSection( int Zelement, double Ein, double outEphoton ) const;

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

FQUALIFIER void 
GUComptonKleinNishina::BuildTable(G4int Z, 
                                  const G4double xmin, 
                                  const G4double xmax,
                                  const G4int nrow,
                                  const G4int ncol,
                                  G4double **p,
                                  )
{
  
  BuildPdfTable(Z,xmin,xmax,nrow,ncol,p); 

  
  fAliasSampler->BuildAliasTable(p);

}

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

FQUALIFIER double 
GUComptonKleinNishina::CalculateDiffCrossSection( int Zelement, 
                                                  double energy0, 
                                                  double energy1 ) const
{
  // based on Geant4 : G4KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  double E0_m = energy0/electron_mass_c2 ;
  double epsilon = energy1/energy0;

  double onecost = (1.- epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);
  double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}

FQUALIFIER void 
GUComptonKleinNishina::BuildPdfTable(G4int Z, 
                                     const G4double xmin, 
                                     const G4double xmax,
                                     const G4int nrow,
                                     const G4int ncol,
                                     G4double **p,
                                     )
{
  // Build the 2-dimensional probability density function (KleinNishina pdf) 
  // in [xmin,xmax] and the alias table (a) and probability (q) with (ncol-1) 
  // equal probable events, each with likelihood 1/(ncol-1)
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 

  //build pdf  
  G4double dx = (xmax - xmin)/nrow;
  G4double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    G4double x = xo + dx*i;

    double ymin = x/(1+2.0*x/electron_mass_c2);
    double dy = (x - ymin)/(n-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;
    for(int j = 0; j < ncol ; ++j) {
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(0,x,y);
      p[i][j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;
    for(int j = 0; i < ncol ; ++j) {
      p[i][j] *= sum;
    }
  }
}

#endif

