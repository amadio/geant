#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "GUPhotoElectronSauterGavrila.h"
#include <iostream>

#include "backend/Backend.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
GUPhotoElectronSauterGavrila::GUPhotoElectronSauterGavrila(Random_t* states, int threadId) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), // fDeltaX(0.1), 
  // fMinY(1.e-8),  fMaxY(1001.), fDeltaY(0.1),
  fMaxZelement(100),       // Elements up to Z=100
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = new GUAliasSampler(fRandomState, fThreadId,
                                     fMaxZelement,
                                     fMinX, fMaxX,
                                     fNrow, fNcol);

  for( int z= 1 ; z < fMaxZelement; ++z)
  {
    //eventually arguments of BuildTable should be replaced by members of *this
    //  and dropped from the function signature. Same for BuildPdfTable
    BuildOneTable(z, fMinX, fMaxX, fNrow, fNcol);
  }
}

VECPHYS_CUDA_HEADER_BOTH 
GUPhotoElectronSauterGavrila::GUPhotoElectronSauterGavrila(Random_t* states, int threadId,
                                             GUAliasSampler* sampler) 
  :
  fRandomState(states), fThreadId(threadId),
  fMinX(1.e-8),  fMaxX(1000.), // fDeltaX(0.1), 
  // fMinY(1.e-8),  fMaxY(1001.), fDeltaY(0.1),  
  fNrow(100), fNcol(100) 
{
  //replace hard coded numbers by default constants
  fAliasSampler = sampler;
}

//need another Ctor with setable parameters

VECPHYS_CUDA_HEADER_BOTH 
GUPhotoElectronSauterGavrila::~GUPhotoElectronSauterGavrila() 
{
  if(fAliasSampler) delete fAliasSampler;
}


VECPHYS_CUDA_HEADER_HOST void 
GUPhotoElectronSauterGavrila::BuildOneTable( int Z, 
                                   const double xmin, 
                                   const double xmax,
                                   const int nrow,
                                   const int ncol)
{
  //for now, the model does not own pdf.  Otherwise, pdf should be the 
  //data member of *this and set its point to the fpdf of fAliasSampler 
  double *pdf = new double [(nrow+1)*ncol];

  BuildLogPdfTable(Z,xmin,xmax,nrow,ncol,pdf); 
  fAliasSampler->BuildAliasTable(Z,nrow,ncol,pdf);

  delete [] pdf;
}

VECPHYS_CUDA_HEADER_HOST void 
GUPhotoElectronSauterGavrila::BuildPdfTable(int Z, 
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
  //  double xo =  xmin + 0.5*dx;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = dx*i;

    const double ymin = -1.0;
    const double dy = 2.0/(ncol-1);
    const double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSectionK(0,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

VECPHYS_CUDA_HEADER_HOST void 
GUPhotoElectronSauterGavrila::BuildLogPdfTable(int Z, 
                                        const double xmin, 
                                        const double xmax,
                                        const int nrow,
                                        const int ncol,
                                        double *p)
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

  //build pdf: logarithmic scale in the input energy bin  

  double logxmin = log(xmin);
  double dx = (log(xmax) - logxmin)/nrow;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    const double ymin = -1.0;
    const double dy = 2./(ncol-1);
    const double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSectionK(Z,x,y);
      p[i*ncol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < ncol ; ++j) {
      p[i*ncol+j] *= sum;
    }
  }
}

// function implementing the angular distribution of photoelectrons

VECPHYS_CUDA_HEADER_BOTH double 
GUPhotoElectronSauterGavrila::CalculateDiffCrossSectionK( int Zelement, 
                                                         double energy, 
                                                         double cosTheta ) const
{
  // based on Geant4 : G4SauterGavrilaAngularDistribution
  // input  : energy   (incomming photon energy)
  //          cosTheta (cons(theta) of photo-electron)
  // output : dsigmaK  (differential cross section, K-shell only) 

  double tau = energy/electron_mass_c2 ;

  double g         = tau + 1.0;
  double invgamma  = 1.0/(tau + 1.0);
  double beta      = sqrt(tau*(tau + 2.0))*invgamma;

  double z  = 1-beta*cosTheta;
  double z2 = z*z;
  double z4 = z2*z2;
  double y  = 1-cosTheta*cosTheta;

  double dsigmaK = (y/z4)*(1+0.5*g*(g-1)*(g-2));

  return dsigmaK;
}

VECPHYS_CUDA_HEADER_BOTH double 
GUPhotoElectronSauterGavrila::CalculateDiffCrossSection( int Zelement, 
                                                         double energy, 
                                                         double cosTheta ) const
{
  // based on Geant4 : G4SauterGavrilaAngularDistribution
  // input  : energy   (incomming photon energy)
  //          cosTheta (cons(theta) of photo-electron)
  // output : dsigma  (differential cross section, K-shell + L-shells) 

  double tau = energy/electron_mass_c2 ;

  double g         = tau + 1.0;
  double invgamma  = 1.0/(tau + 1.0);
  double beta      = sqrt(tau*(tau + 2.0))*invgamma;

  double g2 = g*g;
  double g3 = g2*g;
  double g4 = g2*g2;
  //  double g5 = g2*g3;

  double term = log(g*(1.+beta))/(g*beta);

  double sigmaL2 =    g3 - 5.*g2 + 24.*g - 16. + (g2 + 3*g - 8)*term ;
  double sigmaL3 = 4.*g3 - 6.*g2 +  5.*g +  3. + (g2 - 3*g + 4)*term ;

  double sigma23 = sigmaL2 + sigmaL3;

  double JK = 125./Zelement + 3.5;
  double JL =          1.2;

  double R2 = sigmaL2/sigma23;
  double R3 = sigmaL3/sigma23;

  double PK  = 1 - 1./JK;
  double PL  = (1.-PK)*(1 - 1./JL);
  double PL2 = (1.-PK-PL)*R2;
  double PL3 = (1.-PK-PL)*R3;

  double z  = 1-beta*cosTheta;
  double z2 = z*z;
  double z3 = z2*z;
  double z4 = z2*z2;
  double z5 = z2*z3;
  double y  = 1- cosTheta*cosTheta;

  double dsigmaK = (y/z4)*(1+0.5*g*(g-1)*(g-2))*PK;
  double dsigmaL1 = dsigmaK;

  double coeff= sqrt((g+1)*tau)/pow(g*tau,5.0);

  double dsigmaL2 =  g*(3.*g+1)/(2*z4) - g2*(9*g2+30*g -7)/(8*z3)
                  + g3*(g3+6*g2 +11*g -2)/(4*z2) - g4*(tau*(g+7))/(8*z)
                  + y*((g+1)/z5 - g*(g+1)/z4 + g2*(3*g+1)*(g2-1)/z3);

  dsigmaL2 *= coeff;
  
  double dsigmaL3 = g*(1-3*g)/(2*z4) + g2*(3*g2-1)/z3 
    + g3*(g3-3*g2+2*g+1)/z3 - g4*(g-2)*(g-1)/(2*z)
    +y*( (g+1)/z5 - g*(g+1)*(2*g-1)/z4 + g2*(g2-1)/z3  );

  dsigmaL3 *= coeff;

  double dsigma = PK*dsigmaK + PL*dsigmaL1 +PL2*dsigmaL2 + PL3*dsigmaL3; 
  return dsigma;

}

} // end namespace impl
} // end namespace vecphys

