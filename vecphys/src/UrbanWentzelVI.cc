#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "UrbanWentzelVI.h"
#include <iostream>

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_CUDA_HOST
UrbanWentzelVI::UrbanWentzelVI(Random_t *states, int tid) : EmModelBase<UrbanWentzelVI>(states, tid)
{
  SetLowEnergyLimit(10. * keV);
  SetHighEnergyLimit(100. * MeV);
  Initialization();
}

VECCORE_CUDA_HOST_DEVICE
UrbanWentzelVI::UrbanWentzelVI(Random_t *states, int tid, GUAliasSampler *sampler)
    : EmModelBase<UrbanWentzelVI>(states, tid, sampler)
{
  SetLowEnergyLimit(10. * keV);
  SetHighEnergyLimit(100. * MeV);
}

VECCORE_CUDA_HOST
UrbanWentzelVI::~UrbanWentzelVI() { delete fAliasSampler; }

VECCORE_CUDA_HOST void UrbanWentzelVI::BuildCrossSectionTablePerAtom(int /*Z*/)
{
  ; // dummy for now
}

VECCORE_CUDA_HOST double UrbanWentzelVI::CrossSectionPerAtom(int AtomicNumber, double KineticEnergy)
{
  // G4UrbanMscModel::ComputeCrossSectionPerAtom (Geant4 version 10.2.p1) modified for GeantV
  // (e- MSC for vecphys)

  static const double Tlim = 10.; // in [MeV]
  static const double sigmafactor = twopi * classic_electr_radius * classic_electr_radius;

  static const double epsfactor = 2. * electron_mass_c2 * electron_mass_c2 * Bohr_radius * Bohr_radius / hbarc * hbarc;

  static const double beta2lim =
      Tlim * (Tlim + 2. * electron_mass_c2) / ((Tlim + electron_mass_c2) * (Tlim + electron_mass_c2));
  static const double bg2lim = Tlim * (Tlim + 2. * electron_mass_c2) / (electron_mass_c2 * electron_mass_c2);

  static const double sig0[15] = {0.267, 0.5922, 2.653, 6.235, 11.69, 13.24, 16.12, 23.00,
                                  35.13, 39.95,  50.85, 67.19, 91.15, 104.4, 113.1}; // in [barn]

  static const double Tdat[22] = {0.0001, 0.0002, 0.0004, 0.0007, 0.001, 0.002, 0.004, 0.007, 0.01, 0.02, 0.04, 0.07,
                                  0.1,    0.2,    0.4,    0.7,    1.,    2.,    4.,    7.,    10.,  20.}; // in [MeV]

  static const double epsmin = 1.e-4, epsmax = 1.e10;

  static const double Zdat[15] = {4., 6., 13., 20., 26., 29., 32., 38., 47., 50., 56., 64., 74., 79., 82.};

  // corr. factors for e-/e+ lambda for T <= Tlim
  static const double celectron[15]
                               [22] = {{1.125, 1.072, 1.051, 1.047, 1.047, 1.050, 1.052, 1.054, 1.054, 1.057, 1.062,
                                        1.069, 1.075, 1.090, 1.105, 1.111, 1.112, 1.108, 1.100, 1.093, 1.089, 1.087},
                                       {1.408, 1.246, 1.143, 1.096, 1.077, 1.059, 1.053, 1.051, 1.052, 1.053, 1.058,
                                        1.065, 1.072, 1.087, 1.101, 1.108, 1.109, 1.105, 1.097, 1.090, 1.086, 1.082},
                                       {2.833, 2.268, 1.861, 1.612, 1.486, 1.309, 1.204, 1.156, 1.136, 1.114, 1.106,
                                        1.106, 1.109, 1.119, 1.129, 1.132, 1.131, 1.124, 1.113, 1.104, 1.099, 1.098},
                                       {3.879, 3.016, 2.380, 2.007, 1.818, 1.535, 1.340, 1.236, 1.190, 1.133, 1.107,
                                        1.099, 1.098, 1.103, 1.110, 1.113, 1.112, 1.105, 1.096, 1.089, 1.085, 1.098},
                                       {6.937, 4.330, 2.886, 2.256, 1.987, 1.628, 1.395, 1.265, 1.203, 1.122, 1.080,
                                        1.065, 1.061, 1.063, 1.070, 1.073, 1.073, 1.070, 1.064, 1.059, 1.056, 1.056},
                                       {9.616, 5.708, 3.424, 2.551, 2.204, 1.762, 1.485, 1.330, 1.256, 1.155, 1.099,
                                        1.077, 1.070, 1.068, 1.072, 1.074, 1.074, 1.070, 1.063, 1.059, 1.056, 1.052},
                                       {11.72, 6.364, 3.811, 2.806, 2.401, 1.884, 1.564, 1.386, 1.300, 1.180, 1.112,
                                        1.082, 1.073, 1.066, 1.068, 1.069, 1.068, 1.064, 1.059, 1.054, 1.051, 1.050},
                                       {18.08, 8.601, 4.569, 3.183, 2.662, 2.025, 1.646, 1.439, 1.339, 1.195, 1.108,
                                        1.068, 1.053, 1.040, 1.039, 1.039, 1.039, 1.037, 1.034, 1.031, 1.030, 1.036},
                                       {18.22, 10.48, 5.333, 3.713, 3.115, 2.367, 1.898, 1.631, 1.498, 1.301, 1.171,
                                        1.105, 1.077, 1.048, 1.036, 1.033, 1.031, 1.028, 1.024, 1.022, 1.021, 1.024},
                                       {14.14, 10.65, 5.710, 3.929, 3.266, 2.453, 1.951, 1.669, 1.528, 1.319, 1.178,
                                        1.106, 1.075, 1.040, 1.027, 1.022, 1.020, 1.017, 1.015, 1.013, 1.013, 1.020},
                                       {14.11, 11.73, 6.312, 4.240, 3.478, 2.566, 2.022, 1.720, 1.569, 1.342, 1.186,
                                        1.102, 1.065, 1.022, 1.003, 0.997, 0.995, 0.993, 0.993, 0.993, 0.993, 1.011},
                                       {22.76, 20.01, 8.835, 5.287, 4.144, 2.901, 2.219, 1.855, 1.677, 1.410, 1.224,
                                        1.121, 1.073, 1.014, 0.986, 0.976, 0.974, 0.972, 0.973, 0.974, 0.975, 0.987},
                                       {50.77, 40.85, 14.13, 7.184, 5.284, 3.435, 2.520, 2.059, 1.837, 1.512, 1.283,
                                        1.153, 1.091, 1.010, 0.969, 0.954, 0.950, 0.947, 0.949, 0.952, 0.954, 0.963},
                                       {65.87, 59.06, 15.87, 7.570, 5.567, 3.650, 2.682, 2.182, 1.939, 1.579, 1.325,
                                        1.178, 1.108, 1.014, 0.965, 0.947, 0.941, 0.938, 0.940, 0.944, 0.946, 0.954},
                                       {55.60, 47.34, 15.92, 7.810, 5.755, 3.767, 2.760, 2.239, 1.985, 1.609, 1.343,
                                        1.188, 1.113, 1.013, 0.960, 0.939, 0.933, 0.930, 0.933, 0.936, 0.939, 0.949}};

  // data/corrections for T > Tlim
  static const double hecorr[15] = {120.70, 117.50, 105.00, 92.92, 79.23,  74.510, 68.29, 57.39,
                                    41.97,  36.14,  24.53,  10.21, -7.855, -16.84, -22.30};

  double sigma;

  //  SetParticle(part);

  //  Z23 = G4Pow::GetInstance()->Z23(G4lrint(AtomicNumber));
  double lnZ = math::Log(1.0 * AtomicNumber);
  double Z13 = math::Exp(lnZ / 3.);
  double Z23 = Z13 * Z13;

  // correction if particle .ne. e-/e+
  // compute equivalent kinetic energy
  // lambda depends on p*beta ....

  double eKineticEnergy = KineticEnergy;

  double eTotalEnergy = eKineticEnergy + electron_mass_c2;
  double beta2 = eKineticEnergy * (eTotalEnergy + electron_mass_c2) / (eTotalEnergy * eTotalEnergy);
  double bg2 = eKineticEnergy * (eTotalEnergy + electron_mass_c2) / (electron_mass_c2 * electron_mass_c2);

  double eps = epsfactor * bg2 / Z23;

  if (eps < epsmin)
    sigma = 2. * eps * eps;
  else if (eps < epsmax)
    sigma = math::Log(1. + 2. * eps) - 2. * eps / (1. + 2. * eps);
  else
    sigma = math::Log(2. * eps) - 1. + 1. / eps;

  sigma *= AtomicNumber * AtomicNumber / (beta2 * bg2);

  // interpolate in AtomicNumber and beta2
  double c1, c2, cc1, cc2, corr;

  // get bin number in Z
  int iZ = 14;
  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  while ((iZ >= 0) && (Zdat[iZ] >= AtomicNumber))
    iZ -= 1;
  if (iZ == 14)
    iZ = 13;
  if (iZ == -1)
    iZ = 0;

  double ZZ1 = Zdat[iZ];
  double ZZ2 = Zdat[iZ + 1];
  double ratZ = (AtomicNumber - ZZ1) * (AtomicNumber + ZZ1) / ((ZZ2 - ZZ1) * (ZZ2 + ZZ1));

  if (eKineticEnergy <= Tlim) {
    // get bin number in T (beta2)
    int iT = 21;
    // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
    while ((iT >= 0) && (Tdat[iT] >= eKineticEnergy))
      iT -= 1;
    if (iT == 21)
      iT = 20;
    if (iT == -1)
      iT = 0;

    //  calculate betasquare values
    double T = Tdat[iT], E = T + electron_mass_c2;
    double b2small = T * (E + electron_mass_c2) / (E * E);

    T = Tdat[iT + 1];
    E = T + electron_mass_c2;
    double b2big = T * (E + electron_mass_c2) / (E * E);
    double ratb2 = (beta2 - b2small) / (b2big - b2small);

    c1 = celectron[iZ][iT];
    c2 = celectron[iZ + 1][iT];
    cc1 = c1 + ratZ * (c2 - c1);

    c1 = celectron[iZ][iT + 1];
    c2 = celectron[iZ + 1][iT + 1];
    cc2 = c1 + ratZ * (c2 - c1);

    corr = cc1 + ratb2 * (cc2 - cc1);

    sigma *= sigmafactor / corr;
  } else {
    c1 = bg2lim * sig0[iZ] * (1. + hecorr[iZ] * (beta2 - beta2lim)) / bg2;
    c2 = bg2lim * sig0[iZ + 1] * (1. + hecorr[iZ + 1] * (beta2 - beta2lim)) / bg2;
    if ((AtomicNumber >= ZZ1) && (AtomicNumber <= ZZ2))
      sigma = c1 + ratZ * (c2 - c1);
    else if (AtomicNumber < ZZ1)
      sigma = AtomicNumber * AtomicNumber * c1 / (ZZ1 * ZZ1);
    else if (AtomicNumber > ZZ2)
      sigma = AtomicNumber * AtomicNumber * c2 / (ZZ2 * ZZ2);
    sigma *= barn;
  }
  return sigma;
}

VECCORE_CUDA_HOST void UrbanWentzelVI::Initialization()
{
  if (fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 100, 200);
    BuildAliasTable();
  }
}

VECCORE_CUDA_HOST void UrbanWentzelVI::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function in the input energy randge [fMinX,fMaxX]
  // with an equallogarithmic bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  // output :  p[fNrow][ncol] (probability distribution)
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  const int nrow = fAliasSampler->GetNumEntries();
  const int ncol = fAliasSampler->GetSamplesPerEntry();

  double logxmin = math::Log(fAliasSampler->GetIncomingMin());
  double dx = (math::Log(fAliasSampler->GetIncomingMax()) - logxmin) / nrow;

  for (int i = 0; i <= nrow; ++i) {
    // for each input energy bin
    double x = math::Exp(logxmin + dx * i);

    double ymin = electron_mass_c2;
    double ymax = x - electron_mass_c2;

    double dy = (ymax - ymin) / (ncol - 1);
    double yo = ymin + 0.5 * dy;

    double sum = 0.;

    for (int j = 0; j < ncol; ++j) {
      // for each output energy bin
      double y = yo + dy * j;
      double xsec = CalculateDiffCrossSection(Z, x, y);
      p[i * ncol + j] = xsec;
      sum += xsec;
    }

    // normalization
    const double inv_sum = 1.0 / sum;

    for (int j = 0; j < ncol; ++j) {
      p[i * ncol + j] *= inv_sum;
    }
  }
}

// function implementing the cross section for KleinNishina

VECCORE_CUDA_HOST_DEVICE double UrbanWentzelVI::CalculateDiffCrossSection(int /*Zelement*/, // not used
                                                                          double energy0, double energy1) const
{
  double dsigma = energy1 / energy0;
  return dsigma;
}

VECCORE_CUDA_HOST double UrbanWentzelVI::GetG4CrossSection(int AtomicNumber, double KineticEnergy)
{
  // MSC cross section for the discrete process

  // shut up compiler warnings
  (void) AtomicNumber; (void) KineticEnergy;

  return 0.0;
}

VECCORE_CUDA_HOST_DEVICE void UrbanWentzelVI::SampleByCompositionRejection(int /*Z*/, // not used
                                                                           double energyIn, double &energyOut,
                                                                           double &sinTheta)
{
  energyOut = energyIn;
  sinTheta = 1.0;
}

} // end namespace impl
} // end namespace vecphys
