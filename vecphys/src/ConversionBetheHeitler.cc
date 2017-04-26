#include "ConversionBetheHeitler.h"
#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST ConversionBetheHeitler::ConversionBetheHeitler(Random_t *states, int tid)
    : EmModelBase<ConversionBetheHeitler>(states, tid)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(2. * electron_mass_c2);
  //  Initialization();
}

VECCORE_ATT_HOST_DEVICE ConversionBetheHeitler::ConversionBetheHeitler(Random_t *states, int tid,
                                                                       GUAliasSampler *sampler)
    : EmModelBase<ConversionBetheHeitler>(states, tid, sampler)
{
  fAtomicDependentModel = true;
  SetLowEnergyLimit(2. * electron_mass_c2);
}

VECCORE_ATT_HOST void ConversionBetheHeitler::Initialization()
{
  if (fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 100, 100);
    BuildAliasTable(fAtomicDependentModel);
  }
}

VECCORE_ATT_HOST void ConversionBetheHeitler::BuildCrossSectionTablePerAtom(int /*Z*/)
{
  ; // dummy for now
}

VECCORE_ATT_HOST double ConversionBetheHeitler::GetG4CrossSection(int Z, double GammaEnergy)
{
  // G4BetheHeitlerModel::ComputeCrossSectionPerAtom

  const double GammaEnergyLimit = 1.5 * MeV;

  // Calculates the microscopic cross section in GEANT4 internal units.
  // A parametrized formula from L. Urban is used to estimate
  // the total cross section.
  // It gives a good description of the data from 1.5 MeV to 100 GeV.
  // below 1.5 MeV: sigma=sigma(1.5MeV)*(GammaEnergy-2electronmass)
  //                                   *(GammaEnergy-2electronmass)
  double xSection = 0.0;
  if (Z < 0.9 || GammaEnergy <= 2.0 * electron_mass_c2) {
    return xSection;
  }

  double GammaEnergySave = GammaEnergy;
  if (GammaEnergy < GammaEnergyLimit) {
    GammaEnergy = GammaEnergyLimit;
  }

  double X = math::Log(GammaEnergy / electron_mass_c2), X2 = X * X, X3 = X2 * X, X4 = X3 * X, X5 = X4 * X;

  double F1 = a0 + a1 * X + a2 * X2 + a3 * X3 + a4 * X4 + a5 * X5,
         F2 = b0 + b1 * X + b2 * X2 + b3 * X3 + b4 * X4 + b5 * X5,
         F3 = c0 + c1 * X + c2 * X2 + c3 * X3 + c4 * X4 + c5 * X5;

  xSection = (Z + 1.) * (F1 * Z + F2 * Z * Z + F3) * microbarn;

  if (GammaEnergySave < GammaEnergyLimit) {

    X = (GammaEnergySave - 2. * electron_mass_c2) / (GammaEnergyLimit - 2. * electron_mass_c2);
    xSection *= X * X;
  }

  if (xSection < 0.) {
    xSection = 0.;
  }
  return xSection;
}

VECCORE_ATT_HOST void ConversionBetheHeitler::BuildPdfTable(int Z, double *p)
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

  double logxmin = math::Log(fAliasSampler->GetIncomingMin());
  double dx      = (math::Log(fAliasSampler->GetIncomingMax()) - logxmin) / nrow;

  for (int i = 0; i <= nrow; ++i) {
    // for each input energy bin
    double x = math::Exp(logxmin + dx * i);

    double ymin = electron_mass_c2;
    double ymax = x - electron_mass_c2;

    double dy = (ymax - ymin) / ncol;
    double yo = ymin + 0.5 * dy;

    double sum = 0.;

    for (int j = 0; j < ncol; ++j) {
      // for each output energy bin
      double y        = yo + dy * j;
      double xsec     = CalculateDiffCrossSection(Z, x, y);
      p[i * ncol + j] = xsec;
      sum += xsec;
    }

    // normalization
    sum = 1.0 / sum;

    for (int j = 0; j < ncol; ++j) {
      p[i * ncol + j] *= sum;
    }
  }
}

// function implementing the cross section for KleinNishina
// TODO: need to get electron properties from somewhere

VECCORE_ATT_HOST_DEVICE double ConversionBetheHeitler::CalculateDiffCrossSection(int Zelement, double gammaEnergy,
                                                                                 double electEnergy)
{
  // based on Geant4 : G4BetheHeitlerModel
  // input  : gammaEnergy (incomming photon energy)
  //          electEnergy (converted electron/positron energy)
  // output : dsigma  (differential cross section)

  double epsil = electEnergy / gammaEnergy;

  double epsil0 = electron_mass_c2 / gammaEnergy;
  if (epsil0 > 1.0) {
    return 0;
  }

  // Extract Coulomb factor for this Element
  // F(Z)
  int logZ3 = math::Log(1.0 * int(Zelement + 0.5)) / 3.0;
  double FZ = 8. * logZ3; // 8.*(anElement->GetIonisation()->GetlogZ3());
  if (gammaEnergy > 50. /* *MeV */) {
    FZ += 8. * ComputeCoulombFactor(1.0 * Zelement);
  }

  // delta -> screenvar
  int Z3           = math::Pow(1.0 * int(Zelement + 0.5), 1 / 3.0);
  double screenfac = 136. * epsil0 / Z3; //(anElement->GetIonisation()->GetZ3());
  double screenvar = screenfac / (epsil * (1 - epsil));

  double dsigma = ScreenFunction1(screenvar) * (epsil * epsil + (1. - epsil) * (1. - epsil)) +
                  ScreenFunction2(screenvar) * (2.0 / 3) * epsil * (1.0 - epsil);

  return dsigma;
}

VECCORE_ATT_HOST_DEVICE double ConversionBetheHeitler::ScreenFunction1(double screenVariable) const
{
  // compute the value of the screening function 3*PHI1 - PHI2
  double screenVal;

  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368 * math::Log(screenVariable + 0.952);
  else
    screenVal = 42.392 - screenVariable * (7.796 - 1.961 * screenVariable);

  return screenVal;
}

VECCORE_ATT_HOST_DEVICE double ConversionBetheHeitler::ScreenFunction2(double screenVariable) const
{
  // compute the value of the screening function 1.5*PHI1 - 0.5*PHI2
  double screenVal;

  if (screenVariable > 1.)
    screenVal = 42.24 - 8.368 * math::Log(screenVariable + 0.952);
  else
    screenVal = 41.405 - screenVariable * (5.828 - 0.8945 * screenVariable);

  return screenVal;
}

VECCORE_ATT_HOST_DEVICE
void ConversionBetheHeitler::SampleByCompositionRejection(int elementZ, double GammaEnergy, double &energyOut,
                                                          double &sinTheta)
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

  double epsil;
  double epsil0 = electron_mass_c2 / GammaEnergy;
  if (epsil0 > 1.0) {
    return;
  }

  // do it fast if GammaEnergy < Egsmall
  // select randomly one element constituing the material - input
  const double Egsmall = 2. * MeV;

  if (GammaEnergy < Egsmall) {

    epsil = epsil0 + (0.5 - epsil0) * UniformRandom<double>(fRandomState, fThreadId);
  } else {
    // now comes the case with GammaEnergy >= 2. MeV

    // Extract Coulomb factor for this Element

    double logZ3 = math::Log(1.0 * int(elementZ + 0.5)) / 3.0;
    double FZ    = 8. * logZ3; //(anElement->GetIonisation()->GetlogZ3());
    if (GammaEnergy > 50. * MeV) {
      FZ += 8. * ComputeCoulombFactor(elementZ);
    }

    // limits of the screening variable
    double Z3        = math::Pow(1.0 * int(elementZ + 0.5), 1 / 3.0);
    double screenfac = 136. * epsil0 / Z3; //(anElement->GetIonisation()->GetZ3());
    double screenmax = exp((42.24 - FZ) / 8.368) - 0.952;
    double screenmin = math::Min(4. * screenfac, screenmax);

    // limits of the energy sampling
    double epsil1   = 0.5 - 0.5 * math::Sqrt(1. - screenmin / screenmax);
    double epsilmin = math::Max(epsil0, epsil1), epsilrange = 0.5 - epsilmin;

    //
    // sample the energy rate of the created electron (or positron)
    //
    // double epsil, screenvar, greject ;
    double screenvar, greject;

    double F10    = ScreenFunction1(screenmin) - FZ;
    double F20    = ScreenFunction2(screenmin) - FZ;
    double NormF1 = math::Max(F10 * epsilrange * epsilrange, 0.);
    double NormF2 = math::Max(1.5 * F20, 0.);

    do {
      if (NormF1 / (NormF1 + NormF2) > UniformRandom<double>(fRandomState, fThreadId)) {
        epsil     = 0.5 - epsilrange * math::Pow(UniformRandom<double>(fRandomState, fThreadId), 0.333333);
        screenvar = screenfac / (epsil * (1 - epsil));
        greject   = (ScreenFunction1(screenvar) - FZ) / F10;
      } else {
        epsil     = epsilmin + epsilrange * UniformRandom<double>(fRandomState, fThreadId);
        screenvar = screenfac / (epsil * (1 - epsil));
        greject   = (ScreenFunction2(screenvar) - FZ) / F20;
      }

    } while (greject < UniformRandom<double>(fRandomState, fThreadId));
  } //  end of epsil sampling

  //
  // fixe charges randomly
  //

  double ElectTotEnergy; // PositTotEnergy;
  if (UniformRandom<double>(fRandomState, fThreadId) > 0.5) {
    ElectTotEnergy = (1. - epsil) * GammaEnergy;
    //    PositTotEnergy = epsil*GammaEnergy;
  } else {
    //    PositTotEnergy = (1.-epsil)*GammaEnergy;
    ElectTotEnergy = epsil * GammaEnergy;
  }

  //
  // scattered electron (positron) angles. ( Z - axis along the parent photon)
  //
  //  universal distribution suggested by L. Urban
  // (Geant3 manual (1993) Phys211),
  //  derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  double u;
  // static
  const double aa1 = 0.625;
  const double aa2 = 1.875;
  const double d   = 27.;

  if (9. / (9. + d) > UniformRandom<double>(fRandomState, fThreadId))
    u = -math::Log(UniformRandom<double>(fRandomState, fThreadId) * UniformRandom<double>(fRandomState, fThreadId)) /
        aa1;
  else
    u = -math::Log(UniformRandom<double>(fRandomState, fThreadId) * UniformRandom<double>(fRandomState, fThreadId)) /
        aa2;

  double TetEl = u * electron_mass_c2 / ElectTotEnergy;

  // return energy and sinTheta of the electron -
  // ToDo: store secondaries into a global stack

  energyOut = ElectTotEnergy;
  sinTheta  = math::Sin(TetEl);
}

} // end namespace impl
} // end namespace vecphys
