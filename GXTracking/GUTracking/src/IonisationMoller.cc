#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "IonisationMoller.h"

#include "backend/Backend.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST
IonisationMoller::IonisationMoller(Random_t* states, int tid) 
  : EmModelBase<IonisationMoller>(states,tid)
{
  fDeltaRayThreshold = 1.0*MeV; //temporary: should be set from a cut table 
  SetLowEnergyLimit(0.1*keV);
  BuildAliasTable();
}

VECPHYS_CUDA_HEADER_BOTH 
IonisationMoller::IonisationMoller(Random_t* states, int tid,
                                   GUAliasSampler* sampler) 
  :  EmModelBase<IonisationMoller>(states,tid,sampler)
{
  SetLowEnergyLimit(0.1*keV);
}

VECPHYS_CUDA_HEADER_HOST void 
IonisationMoller::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECPHYS_CUDA_HEADER_HOST void 
IonisationMoller::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (MollerBhabha pdf) in the
  // input energy randge [fMinX,fMaxX] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  // output :  p[fNrow][fNcol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  double logxmin = log(fMinX);
  double dx = (log(fMaxX) - logxmin)/fNrow;
  //  double xo =  fMinX + 0.5*dx;

  for(int i = 0; i <= fNrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    //e-e- (Moller) only for now
    double ymin = fDeltaRayThreshold; //minimum delta-ray energy
    double dy = (x/2.0 - ymin)/fNcol; //maximum x/2.0
    double yo = ymin + 0.5*dy;
  
    double sum = 0.;

    for(int j = 0; j < fNcol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
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

// function implementing the cross section for MollerBhabha

VECPHYS_CUDA_HEADER_BOTH double 
IonisationMoller::CalculateDiffCrossSection(int Zelement,
                                            double kineticEnergy, 
					    double deltaRayEnergy) const
{
  // based on Geant3 : Simulation of the delta-ray production (PHY331-1)
  // input  : kineticEnergy (incomming photon energy)
  //          deltaRayEnergy (scattered photon energy)
  // output : dcross (differential cross section) 

  double dcross = 0.0;

  double tau = kineticEnergy/electron_mass_c2;
  double gam = tau + 1.0;
  double gamma2 = gam*gam;

  double epsil = deltaRayEnergy/kineticEnergy;

  //Moller (e-e-) scattering only
  //Bhabha (e+e-) scattering not implemented for now

  double fgam = (2*gam-1.0)/gamma2;
  double x = 1/epsil;
  double y = 1/(1.0-epsil);

  dcross = 1-fgam + x*(x-fgam) + y*(y-fgam);

  return dcross;
}

VECPHYS_CUDA_HEADER_BOTH double 
IonisationMoller::GetG4CrossSection(double kineticEnergy, 
                                    const int Z)
{
  G4double cross = 0.0;

  //temporary - set by material
  G4double cutEnergy = fDeltaRayThreshold;
  G4double maxEnergy = 1.0*TeV;

  G4double tmax =  0.5*kineticEnergy;
  tmax = Min(maxEnergy, tmax);

  if(cutEnergy < tmax) {
    G4double xmin  = cutEnergy/kineticEnergy;
    G4double xmax  = tmax/kineticEnergy;
    G4double tau   = kineticEnergy/electron_mass_c2;
    G4double gam   = tau + 1.0;
    G4double gamma2= gam*gam;
    G4double beta2 = tau*(tau + 2)/gamma2;

    //Moller (e-e-) scattering

    G4double gg = (2.0*gam - 1.0)/gamma2;
    cross = ((xmax - xmin)*(1.0 - gg + 1.0/(xmin*xmax)
          + 1.0/((1.0-xmin)*(1.0 - xmax)))
          - gg*G4Log( xmax*(1.0 - xmin)/(xmin*(1.0 - xmax)) ) ) / beta2;

  }
  cross *= Z*twopi_mc2_rcl2/kineticEnergy;

  return cross;
}

VECPHYS_CUDA_HEADER_BOTH void 
IonisationMoller::SampleByCompositionRejection(int     Z, //not used
                                               double  kineticEnergy,
                                               double& deltaKinEnergy,
                                               double& sinTheta)
{
  //temporary - set by material
  G4double cutEnergy = fDeltaRayThreshold;
  G4double maxEnergy = 1.0*TeV;

  //based on G4MollerBhabhaModel::SampleSecondaries

  G4double tmin = cutEnergy;  
  G4double tmax = 0.5*kineticEnergy; 

  if(maxEnergy < tmax) { tmax = maxEnergy; }

  if(tmin >= tmax) { return; }

  G4double energy = kineticEnergy + electron_mass_c2;
  G4double xmin   = tmin/kineticEnergy;
  G4double xmax   = tmax/kineticEnergy;
  G4double gam    = energy/electron_mass_c2;
  G4double gamma2 = gam*gam;
  //  G4double beta2  = 1.0 - 1.0/gamma2;
  G4double x, z, q, grej;

  //Moller (e-e-) scattering
  G4double gg = (2.0*gam - 1.0)/gamma2;
  G4double y = 1.0 - xmax;
  grej = 1.0 - gg*xmax + xmax*xmax*(1.0 - gg + (1.0 - gg*y)/(y*y));

  do {
    q = UniformRandom<kScalar>(fRandomState,fThreadId);
    x = xmin*xmax/(xmin*(1.0 - q) + xmax*q);
    y = 1.0 - x;
    z = 1.0 - gg*x + x*x*(1.0 - gg + (1.0 - gg*y)/(y*y));
  } while(grej * UniformRandom<kScalar>(fRandomState,fThreadId) > z);

  deltaKinEnergy = x * kineticEnergy;

  G4double totalMomentum = sqrt(kineticEnergy*(kineticEnergy + electron_mass_c2));

  G4double deltaMomentum =
    sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0*electron_mass_c2));
  G4double cost = deltaKinEnergy * (energy + electron_mass_c2) /
    (deltaMomentum * totalMomentum );
  if(cost > 1.0) { cost = 1.0; }
  G4double sint2 = (1.0 - cost)*(1.0 + cost);

  sinTheta = (sint2 < 0.0) ? 0.0 : sqrt(sint2);
}

} // end namespace impl
} // end namespace vecphys
