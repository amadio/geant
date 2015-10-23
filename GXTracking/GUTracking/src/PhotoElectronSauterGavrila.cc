#include <iostream>

#include "backend/Backend.h"
// #include "GUAuxFunctions.h"
#include "GUG4TypeDef.h"

#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "PhotoElectronSauterGavrila.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
PhotoElectronSauterGavrila::PhotoElectronSauterGavrila(Random_t* states, int tid) 
  : EmModelBase<PhotoElectronSauterGavrila>(states,tid)
{
  SetLowEnergyLimit(10.*keV);
  BuildAliasTable();
}

VECPHYS_CUDA_HEADER_BOTH 
PhotoElectronSauterGavrila::PhotoElectronSauterGavrila(Random_t* states, int tid,
                                                       GUAliasSampler* sampler) 
  : EmModelBase<PhotoElectronSauterGavrila>(states,tid,sampler)
{
  SetLowEnergyLimit(10.*keV);
}

VECPHYS_CUDA_HEADER_HOST void 
PhotoElectronSauterGavrila::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECPHYS_CUDA_HEADER_HOST void 
PhotoElectronSauterGavrila::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (KleinNishina pdf) in the 
  // input energy randge [fMinX,fMaxX] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number)
  // output :  p[fNrow][fNcol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  double logxmin = log(fMinX);
  double dx = (log(fMaxX) - logxmin)/fNrow;

  for(int i = 0; i <= fNrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    const double ymin = -1.0;
    const double dy = 2./(fNcol-1); 
    const double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < fNcol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSectionK(Z,x,y);

      p[i*fNcol+j] = xsec;
      sum += xsec;
    }

    //normalization
    sum = 1.0/sum;

    for(int j = 0; j < fNcol ; ++j) {
      p[i*fNcol+j] *= sum;
    }
  }
}

// function implementing the angular distribution of photoelectrons

VECPHYS_CUDA_HEADER_BOTH double 
PhotoElectronSauterGavrila::CalculateDiffCrossSectionK(int Zelement, 
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
PhotoElectronSauterGavrila::CalculateDiffCrossSection(int Zelement, 
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

VECPHYS_CUDA_HEADER_BOTH double
PhotoElectronSauterGavrila::GetG4CrossSection(double  gammaEnergy, 
                                              const int Z)
{
  double xSection = 0.;
  //dummy for now
  return xSection;
}

VECPHYS_CUDA_HEADER_BOTH void 
PhotoElectronSauterGavrila::SampleByCompositionRejection(int    Z, //not used
                                                         double energyIn,
                                                         double& energyOut,
                                                         double& cost)
{
  //use the scalar implementation which is equivalent to Geant4
  energyOut =  GetPhotoElectronEnergy<kScalar>(energyIn,Z);

  //sample angular direction according to SauterGavrilaAngularDistribution

  G4double tau = energyIn/electron_mass_c2;
  //  static 
  const G4double taulimit = 50.0;

  if (tau > taulimit) {
    cost  = 1.0; //set to the primary direction
  } else {
    // algorithm according Penelope 2008 manual and 
    // F.Sauter Ann. Physik 9, 217(1931); 11, 454(1931). 

    G4double gamma = tau + 1;
    G4double beta  = std::sqrt(tau*(tau + 2))/gamma;
    G4double A     = (1 - beta)/beta;
    G4double Ap2   = A + 2;
    G4double B     = 0.5*beta*gamma*(gamma - 1)*(gamma - 2);
    G4double grej  = 2*(1 + A*B)/A;
    G4double z, g;
    do { 
      G4double q = UniformRandom<kScalar>(fRandomState,fThreadId);
      z = 2*A*(2*q + Ap2*std::sqrt(q))/(Ap2*Ap2 - 4*q); 
      g = (2 - z)*(1.0/(A + z) + B);

    } while(g < UniformRandom<kScalar>(fRandomState,fThreadId)*grej);
 
    cost = 1 - z;
  }
}

} // end namespace impl
} // end namespace vecphys
