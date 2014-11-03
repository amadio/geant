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


FQUALIFIER 
GUComptonKleinNishina::GUComptonKleinNishina() 
{
}

FQUALIFIER void 
GUComptonKleinNishina::BuildTable( int Z, 
                                   const double xmin, 
                                   const double xmax,
                                   const int nrow,
                                   const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should
  //be the data member of *this

  double *pdf = (double*) malloc(nrow*ncol*sizeof(double));

  BuildPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTables(nrow,ncol,pdf);

  free(pdf);
}

FQUALIFIER void 
GUComptonKleinNishina::BuildPdfTable(int Z, 
                                     const double xmin, 
                                     const double xmax,
                                     const int nrow,
                                     const int ncol,
                                     double *p
                                     )
{
  // Build the probability density function (KleinNishina pdf) 
  // in the energy randge [xmin,xmax] with an equal bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  //           xmin (miminum energy)
  //           xmax (maxinum energy)
  //           nrow (number of input energy bins)
  //           ncol (number of output energy bins)
  //
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  //build pdf  
  double dx = (xmax - xmin)/nrow;
  double xo =  xmin + 0.5*dx;

  for(int i = 0; i < nrow ; ++i) {
    //for each input energy bin
    double x = xo + dx*i;

    double ymin = x/(1+2.0*x/electron_mass_c2);
    double dy = (x - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;
    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(0,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;
    for(int j = 0; i < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

FQUALIFIER double 
GUComptonKleinNishina::CalculateDiffCrossSection( int Zelement, 
                                                  double energy0, 
                                                  double energy1 ) const
{
  // based on Geant4 : KleinNishinaCompton
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

#endif

