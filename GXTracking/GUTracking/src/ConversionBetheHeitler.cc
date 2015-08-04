#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "ConversionBetheHeitler.h"

#include "backend/Backend.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECPHYS_CUDA_HEADER_HOST ConversionBetheHeitler::
ConversionBetheHeitler(Random_t* states, int tid) 
  :  EmModelBase<ConversionBetheHeitler>(states,tid)
{
  SetLowEnergyLimit(2.*electron_mass_c2);
  BuildAliasTable();
}

VECPHYS_CUDA_HEADER_BOTH ConversionBetheHeitler::
ConversionBetheHeitler(Random_t* states, int tid,
                         GUAliasSampler* sampler) 
  : EmModelBase<ConversionBetheHeitler>(states,tid,sampler)
{
  SetLowEnergyLimit(2.*electron_mass_c2);
}

VECPHYS_CUDA_HEADER_HOST void 
ConversionBetheHeitler::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECPHYS_CUDA_HEADER_HOST void 
ConversionBetheHeitler::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (BetheHeitler pdf) in the
  // input energy randge [xmin,xmax] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number) 
  // output :  p[nrow][ncol] (probability distribution) 
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  double logxmin = log(fMinX);

  double dx = (log(fMaxX) - logxmin)/fNrow;

  for(int i = 0; i <= fNrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    double ymin = electron_mass_c2;
    double ymax = x - electron_mass_c2;

    double dy = (ymax - ymin)/(fNcol-1);
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

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECPHYS_CUDA_HEADER_BOTH double 
ConversionBetheHeitler::CalculateDiffCrossSection(int Zelement, 
                                                  double gammaEnergy, 
                                                  double electEnergy)
{ 
  // based on Geant4 : G4BetheHeitlerModel
  // input  : gammaEnergy (incomming photon energy)
  //          electEnergy (converted electron/positron energy)
  // output : dsigma  (differential cross section) 

  double epsil  = electEnergy/gammaEnergy;

  double epsil0 = electron_mass_c2/gammaEnergy ;
  if(epsil0 > 1.0) { return 0; }

  // Extract Coulomb factor for this Element
  //F(Z)
  int logZ3 = log(1.0*int(Zelement + 0.5))/3.0; 
  double FZ = 8.*logZ3 ; //8.*(anElement->GetIonisation()->GetlogZ3());
  if (gammaEnergy > 50. /* *MeV */) { 
    FZ += 8.*ComputeCoulombFactor(1.0*Zelement); 
  }
  
  //delta -> screenvar
  int Z3 = pow(1.0*int(Zelement + 0.5),1/3.0);
  double screenfac = 136.*epsil0/Z3; //(anElement->GetIonisation()->GetZ3());
  double screenvar = screenfac/(epsil*(1-epsil));
  
  double dsigma = ScreenFunction1(screenvar)*(epsil*epsil+(1.-epsil)*(1.-epsil))  
                + ScreenFunction2(screenvar)*(2.0/3)*epsil*(1.0-epsil);
  
  return dsigma;
}

	 /*
VECPHYS_CUDA_HEADER_BOTH double 
ConversionBetheHeitler::ComputeCoulombFactor(double fZeff) const
{
  // Compute Coulomb correction factor (Phys Rev. D50 3-1 (1994) page 1254)

  const double k1 = 0.0083 , k2 = 0.20206 ,k3 = 0.0020 , k4 = 0.0369 ;
  const double fine_structure_const = (1.0/137); //check unit

  double az1 = fine_structure_const*fZeff;
  double az2 = az1 * az1;
  double az4 = az2 * az2;

  double fCoulomb = (k1*az4 + k2 + 1./(1.+az2))*az2 - (k3*az4 + k4)*az4;
  return fCoulomb;
}
	 */

VECPHYS_CUDA_HEADER_BOTH double 
ConversionBetheHeitler::ScreenFunction1(double screenVariable) const
{
  // compute the value of the screening function 3*PHI1 - PHI2
  double screenVal;
  
  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368*log(screenVariable+0.952);
  else
    screenVal = 42.392 - screenVariable*(7.796 - 1.961*screenVariable);
  
  return screenVal;
}

VECPHYS_CUDA_HEADER_BOTH double 
ConversionBetheHeitler::ScreenFunction2(double screenVariable) const
{
  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
  double screenVal;
  
  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368*log(screenVariable+0.952);
  else
    screenVal = 41.405 - screenVariable*(5.828 - 0.8945*screenVariable);

  return screenVal;
}

VECPHYS_CUDA_HEADER_BOTH 
void ConversionBetheHeitler::
SampleByCompositionRejection(int     elementZ,
                             double  GammaEnergy,
                             double& energyOut,
                             double& sinTheta)
{
// G4BetheHeitlerModel::SampleSecondaries
//
// The secondaries e+e- energies are sampled using the Bethe - Heitler
// cross sections with Coulomb correction.
// A modified version of the random number techniques of Butcher & Messel
// is used (Nuc Phys 20(1960),15).
//
// GEANT4 internal units.
//
// Note 1 : Effects due to the breakdown of the Born approximation at
//          low energy are ignored.
// Note 2 : The differential cross section implicitly takes account of 
//          pair creation in both nuclear and atomic electron fields.
//          However triplet prodution is not generated.

  G4double epsil ;
  G4double epsil0 = electron_mass_c2/GammaEnergy ;
  if(epsil0 > 1.0) { return; }

  // do it fast if GammaEnergy < Egsmall
  // select randomly one element constituing the material - input
  const double Egsmall=2.*MeV;

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5-epsil0)*UniformRandom<kScalar>(fRandomState,fThreadId);

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element
    int logZ3 = log(1.0*int(elementZ + 0.5))/3.0; 
    G4double FZ = 8.*logZ3; //(anElement->GetIonisation()->GetlogZ3());
    if (GammaEnergy > 50.*MeV) { FZ += 8.*ComputeCoulombFactor(elementZ); }

    // limits of the screening variable
    int Z3 = pow(1.0*int(elementZ + 0.5),1/3.0);
    G4double screenfac = 136.*epsil0/Z3;//(anElement->GetIonisation()->GetZ3());
    G4double screenmax = exp ((42.24 - FZ)/8.368) - 0.952 ;
    G4double screenmin = fmin(4.*screenfac,screenmax);

    // limits of the energy sampling
    G4double epsil1 = 0.5 - 0.5*sqrt(1. - screenmin/screenmax) ;
    G4double epsilmin = fmax(epsil0,epsil1) , epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    //G4double epsil, screenvar, greject ;
    G4double  screenvar, greject ;

    G4double F10 = ScreenFunction1(screenmin) - FZ;
    G4double F20 = ScreenFunction2(screenmin) - FZ;
    G4double NormF1 = fmax(F10*epsilrange*epsilrange,0.); 
    G4double NormF2 = fmax(1.5*F20,0.);

    do {
      if ( NormF1/(NormF1+NormF2) > UniformRandom<kScalar>(fRandomState,fThreadId) ) {
        epsil = 0.5 - epsilrange*pow(UniformRandom<kScalar>(fRandomState,fThreadId), 0.333333);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction1(screenvar) - FZ)/F10;
      } else { 
        epsil = epsilmin + epsilrange*UniformRandom<kScalar>(fRandomState,fThreadId);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction2(screenvar) - FZ)/F20;
      }

    } while( greject < UniformRandom<kScalar>(fRandomState,fThreadId));
  }   //  end of epsil sampling
   
  //
  // fixe charges randomly
  //

  G4double ElectTotEnergy;// PositTotEnergy;
  if ( UniformRandom<kScalar>(fRandomState,fThreadId) > 0.5) {
    ElectTotEnergy = (1.-epsil)*GammaEnergy;
    //    PositTotEnergy = epsil*GammaEnergy;
     
  } else {
    //    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil*GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban 
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double u;
  //static
  const G4double aa1 = 0.625; 
  const G4double aa2 = 1.875;
  const G4double d = 27. ;

  if (9./(9.+d) >UniformRandom<kScalar>(fRandomState,fThreadId)) 
   u= - G4Log(UniformRandom<kScalar>(fRandomState,fThreadId)*UniformRandom<kScalar>(fRandomState,fThreadId))/aa1;
  else                            
    u= - G4Log(UniformRandom<kScalar>(fRandomState,fThreadId)*UniformRandom<kScalar>(fRandomState,fThreadId))/aa2;

  G4double TetEl = u*electron_mass_c2/ElectTotEnergy;

  /*
  G4double TetPo = u*electron_mass_c2/PositTotEnergy;
  G4double Phi  = twopi *UniformRandom<kScalar>(fRandomState,fThreadId) ;
  G4double dxEl= sin(TetEl)*cos(Phi),dyEl= sin(TetEl)*sin(Phi),dzEl=cos(TetEl);
  G4double dxPo=-sin(TetPo)*cos(Phi),dyPo=-sin(TetPo)*sin(Phi),dzPo=cos(TetPo);
  */

  //return energy and sinTheta of the electron - 
  //ToDo: store secondaries into a global stack
  energyOut = ElectTotEnergy;
  sinTheta = sin(TetEl);

}

} // end namespace impl
} // end namespace vecphys