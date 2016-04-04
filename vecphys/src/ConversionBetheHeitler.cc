#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "ConversionBetheHeitler.h"

#include "base/Global.h"
#include "GUG4TypeDef.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST ConversionBetheHeitler::
ConversionBetheHeitler(Random_t* states, int tid)
  : EmModelBase<ConversionBetheHeitler>(states,tid)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(2.*electron_mass_c2);
  Initialization();
}

VECCORE_CUDA_HOST_DEVICE ConversionBetheHeitler::
ConversionBetheHeitler(Random_t* states, int tid,
                         GUAliasSampler* sampler)
  : EmModelBase<ConversionBetheHeitler>(states,tid,sampler)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(2.*electron_mass_c2);
}

VECCORE_CUDA_HOST void
ConversionBetheHeitler::Initialization()
{
  if(fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId,
				       fLowEnergyLimit, fHighEnergyLimit,
                                       100, 100);
    BuildAliasTable(fAtomicDependentModel);
  }
}

VECCORE_CUDA_HOST void
ConversionBetheHeitler::BuildCrossSectionTablePerAtom(int Z)
{
  ; //dummy for now
}

VECCORE_CUDA_HOST void
ConversionBetheHeitler::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (BetheHeitler pdf) in the
  // input energy randge [xmin,xmax] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number)
  // output :  p[nrow][ncol] (probability distribution)
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  const int nrow = fAliasSampler->GetNumEntries();
  const int ncol = fAliasSampler->GetSamplesPerEntry();

  double logxmin = log(fAliasSampler->GetIncomingMin());
  double dx = (log(fAliasSampler->GetIncomingMax()) - logxmin)/nrow;

  for(int i = 0; i <= nrow ; ++i) {
    //for each input energy bin
    double x = exp(logxmin + dx*i);

    double ymin = electron_mass_c2;
    double ymax = x - electron_mass_c2;

    double dy = (ymax - ymin)/(ncol-1);
    double yo = ymin + 0.5*dy;

    double sum = 0.;

    for(int j = 0; j < ncol ; ++j) {
      //for each output energy bin
      double y = yo + dy*j;
      double xsec = CalculateDiffCrossSection(Z,x,y);
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

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECCORE_CUDA_HOST_DEVICE double
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

VECCORE_CUDA_HOST_DEVICE double
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

VECCORE_CUDA_HOST_DEVICE double
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

VECCORE_CUDA_HOST_DEVICE
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

    epsil = epsil0 + (0.5-epsil0)*UniformRandom<backend::Scalar>(fRandomState,fThreadId);

  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element

    G4double logZ3 = log(1.0*int(elementZ + 0.5))/3.0;
    G4double FZ = 8.*logZ3; //(anElement->GetIonisation()->GetlogZ3());
    if (GammaEnergy > 50.*MeV) { FZ += 8.*ComputeCoulombFactor(elementZ); }

    // limits of the screening variable
    G4double Z3 = pow(1.0*int(elementZ + 0.5),1/3.0);
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
      if ( NormF1/(NormF1+NormF2) > UniformRandom<backend::Scalar>(fRandomState,fThreadId) ) {
	epsil = 0.5 - epsilrange*pow(UniformRandom<backend::Scalar>(fRandomState,fThreadId), 0.333333);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction1(screenvar) - FZ)/F10;
      } else {
        epsil = epsilmin + epsilrange*UniformRandom<backend::Scalar>(fRandomState,fThreadId);
        screenvar = screenfac/(epsil*(1-epsil));
        greject = (ScreenFunction2(screenvar) - FZ)/F20;
      }

    } while( greject < UniformRandom<backend::Scalar>(fRandomState,fThreadId));
  }   //  end of epsil sampling

  //
  // fixe charges randomly
  //

  G4double ElectTotEnergy;// PositTotEnergy;
  if ( UniformRandom<backend::Scalar>(fRandomState,fThreadId) > 0.5) {
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

  if (9./(9.+d) >UniformRandom<backend::Scalar>(fRandomState,fThreadId))
   u= - G4Log(UniformRandom<backend::Scalar>(fRandomState,fThreadId)*UniformRandom<backend::Scalar>(fRandomState,fThreadId))/aa1;
  else
    u= - G4Log(UniformRandom<backend::Scalar>(fRandomState,fThreadId)*UniformRandom<backend::Scalar>(fRandomState,fThreadId))/aa2;

  G4double TetEl = u*electron_mass_c2/ElectTotEnergy;

  /*
  G4double TetPo = u*electron_mass_c2/PositTotEnergy;
  G4double Phi  = twopi *UniformRandom<backend::Scalar>(fRandomState,fThreadId) ;
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
