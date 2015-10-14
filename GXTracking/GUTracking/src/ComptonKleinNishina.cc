#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "ComptonKleinNishina.h"
#include <iostream>

#include "backend/Backend.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
ComptonKleinNishina::ComptonKleinNishina(Random_t* states, int tid) 
  : EmModelBase<ComptonKleinNishina>(states,tid)
{
  SetLowEnergyLimit(10.*keV);
  SetNcol(400);

  BuildAliasTable();
}

VECPHYS_CUDA_HEADER_BOTH 
ComptonKleinNishina::ComptonKleinNishina(Random_t* states, int tid,
                                         GUAliasSampler* sampler) 
  : EmModelBase<ComptonKleinNishina>(states,tid,sampler)
{
  SetLowEnergyLimit(10.*keV);
}

VECPHYS_CUDA_HEADER_HOST void 
ComptonKleinNishina::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECPHYS_CUDA_HEADER_HOST void 
ComptonKleinNishina::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (KleinNishina pdf) in the
  // input energy randge [fMinX,fMaxX] with an equallogarithmic  bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  // output :  p[fNrow][fNcol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  double logxmin = log(fMinX);
  double logxmax = log(fMaxX);
  double dx = (logxmax - logxmin)/fNrow;

  int    nintegral= 5*int(log(logxmax)); //temporary
  double average = 0;
  double normal = 0;
  double fxsec = 0;

  for(int i = 0; i <= fNrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    double ymin = x/(1+2.0*x*inv_electron_mass_c2);
    double dy = (x - ymin)/fNcol;
    //    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < fNcol ; ++j) {
      //for each output energy bin
      //      double y = yo + dy*j;
      double y = 0;
      average = 0;
      normal = 0;

      //cross section weighted bin position
      for(int k = 0; k < nintegral ; ++k) {
        y = ymin + dy*(j+(0.5+k)/nintegral);
        fxsec = CalculateDiffCrossSection(Z,x,y);
        average += y*fxsec; 
        normal  += fxsec;
      }
      y = average/normal;

      double xsec = CalculateDiffCrossSection(Z,x,y);
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

// function implementing the cross section for KleinNishina

VECPHYS_CUDA_HEADER_BOTH double 
ComptonKleinNishina::CalculateDiffCrossSection(int Zelement, 
                                               double energy0, 
                                               double energy1 ) const
{
  // based on Geant4 : KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section) 

  double E0_m = energy0*inv_electron_mass_c2 ;
  double epsilon = energy1/energy0;

  double onecost = (1.- epsilon)/(epsilon*E0_m);
  double sint2   = onecost*(2.-onecost);
  double greject = 1. - epsilon*sint2/(1.+ epsilon*epsilon);
  double dsigma = (epsilon + 1./epsilon)*greject;

  return dsigma;
}

VECPHYS_CUDA_HEADER_BOTH double
ComptonKleinNishina::GetG4CrossSection(double  gammaEnergy, 
                                       const int Z)
{
  //G4KleinNishinaModel::ComputeCrossSectionPerAtom - Genat4 10.1.p2
  double xSection = 0.;

  const G4double dT0 = keV;
  const G4double a = 20.0 , b = 230.0 , c = 440.0;

  /*  
  G4double p1Z = Z*(d1 + e1*Z + f1*Z*Z);
  G4double p2Z = Z*(d2 + e2*Z + f2*Z*Z);
  G4double p3Z = Z*(d3 + e3*Z + f3*Z*Z); 
  G4double p4Z = Z*(d4 + e4*Z + f4*Z*Z);
  */

  G4double p1Z = Z*( 2.7965e-1 +  1.9756e-5*Z + -3.9178e-7*Z*Z)*barn;
  G4double p2Z = Z*(-1.8300e-1 + -1.0205e-2*Z +  6.8241e-5*Z*Z)*barn;
  G4double p3Z = Z*( 6.7527    + -7.3913e-2*Z +  6.0480e-5*Z*Z)*barn;
  G4double p4Z = Z*(-1.9798e+1 +  2.7079e-2*Z +  3.0274e-4*Z*Z)*barn;

  G4double T0  = 15.0*keV; 
  if (Z < 1.5) { T0 = 40.0*keV; } 

  G4double X   = Max(gammaEnergy, T0) / electron_mass_c2;
  xSection = p1Z*G4Log(1.+2.*X)/X
               + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    
  //  modification for low energy. (special case for Hydrogen)
  if (gammaEnergy < T0) {
    X = (T0+dT0) / electron_mass_c2 ;
    G4double sigma = p1Z*G4Log(1.+2*X)/X
                    + (p2Z + p3Z*X + p4Z*X*X)/(1. + a*X + b*X*X + c*X*X*X);
    G4double   c1 = -T0*(sigma-xSection)/(xSection*dT0);             
    G4double   c2 = 0.150; 
    if (Z > 1.5) { c2 = 0.375-0.0556*G4Log(1.0*Z); }
    G4double    y = G4Log(gammaEnergy/T0);
    xSection *= G4Exp(-y*(c1+c2*y));          
  }
  if(xSection < 0.0) { xSection = 0.0; }
  return xSection;
}

VECPHYS_CUDA_HEADER_BOTH void 
ComptonKleinNishina::SampleByCompositionRejection(int    Z, //not used
                                                  double energyIn,
                                                  double& energyOut,
                                                  double& sinTheta)
{
  double epsilon, epsilonsq, onecost, sint2, greject ;
  
  double E0_m = energyIn*inv_electron_mass_c2 ;
  double eps0       = 1./(1. + 2.*E0_m);
  double epsilon0sq = eps0*eps0;
  double alpha1     = - log(eps0);
  double alpha2     = 0.5*(1.- epsilon0sq);
  
  do {
    if(alpha1/(alpha1+alpha2) > UniformRandom<kScalar>(fRandomState,fThreadId))
    {
      epsilon   = exp(-alpha1*UniformRandom<kScalar>(fRandomState,fThreadId));
      epsilonsq = epsilon*epsilon; 
    } 
    else {
      epsilonsq = epsilon0sq+(1.- epsilon0sq)*UniformRandom<kScalar>(fRandomState,fThreadId);
      epsilon   = sqrt(epsilonsq);
    }
    
    onecost = (1.- epsilon)/(epsilon*E0_m);
    sint2   = onecost*(2.-onecost);
    greject = 1. - epsilon*sint2/(1.+ epsilonsq);
    
  } while (greject < UniformRandom<kScalar>(fRandomState,fThreadId));
  
  energyOut = epsilon*energyIn;
  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

} // end namespace impl
} // end namespace vecphys
