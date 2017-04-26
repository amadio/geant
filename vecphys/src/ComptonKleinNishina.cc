#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "ComptonKleinNishina.h"
#include <iostream>

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
ComptonKleinNishina::ComptonKleinNishina(Random_t *states, int tid) : EmModelBase<ComptonKleinNishina>(states, tid)
{
  SetLowEnergyLimit(10. * keV);
  //  Initialization();
}

VECCORE_ATT_HOST_DEVICE
ComptonKleinNishina::ComptonKleinNishina(Random_t *states, int tid, GUAliasSampler *sampler)
    : EmModelBase<ComptonKleinNishina>(states, tid, sampler)
{
  SetLowEnergyLimit(10. * keV);
}

VECCORE_ATT_HOST
ComptonKleinNishina::~ComptonKleinNishina() { delete fAliasSampler; }

VECCORE_ATT_HOST void ComptonKleinNishina::BuildCrossSectionTablePerAtom(int /*Z*/)
{
  ; // dummy for now
}

VECCORE_ATT_HOST void ComptonKleinNishina::Initialization()
{
  if (fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 100, 200);
    BuildAliasTable();
  }
}

VECCORE_ATT_HOST void ComptonKleinNishina::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (KleinNishina pdf) in the
  // input energy randge [fMinX,fMaxX] with an equallogarithmic  bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
  // output :  p[fNrow][ncol] (probability distribution)
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  const int nrow = fAliasSampler->GetNumEntries();
  const int ncol = fAliasSampler->GetSamplesPerEntry();

  double logxmin = math::Log(fAliasSampler->GetIncomingMin());
  double logxmax = math::Log(fAliasSampler->GetIncomingMax());
  double dx = (logxmax - logxmin) / nrow;

  int nintegral = 5 * int(math::Log(logxmax)); // temporary

  // use the average xsec within the bin instead of xsec at the mid-point
  for (int i = 0; i <= nrow; ++i) {

    double x = math::Exp(logxmin + dx * i);

    double ymin = x / (1 + 2.0 * x * inv_electron_mass_c2);
    double dy = (x - ymin) / ncol;

    double sum = 0.;
    for (int j = 0; j < ncol; ++j) {
      // for each input energy bin
      double normal = 0; // rename to sumN    // Sum N (pdf)

      // cross section weighted bin position
      for (int k = 0; k < nintegral; ++k) {

        double y = ymin + dy * (j + (0.5 + k) / nintegral);
        double fxsec = CalculateDiffCrossSection(Z, x, y);
        normal += fxsec;
      }

      double xsec = normal / nintegral;
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

VECCORE_ATT_HOST_DEVICE double ComptonKleinNishina::CalculateDiffCrossSection(int /*Zelement*/, // not used
                                                                               double energy0, double energy1) const
{
  // based on Geant4 : KleinNishinaCompton
  // input  : energy0 (incomming photon energy)
  //          energy1 (scattered photon energy)
  // output : dsigma  (differential cross section)

  double E0_m = energy0 * inv_electron_mass_c2;
  double epsilon = energy1 / energy0;

  double onecost = (1. - epsilon) / (epsilon * E0_m);
  double sint2 = onecost * (2. - onecost);
  double greject = 1. - epsilon * sint2 / (1. + epsilon * epsilon);
  double dsigma = (epsilon + 1. / epsilon) * greject;

  return dsigma;
}

VECCORE_ATT_HOST double ComptonKleinNishina::GetG4CrossSection(int Z, double gammaEnergy)
{
  // G4KleinNishinaModel::ComputeCrossSectionPerAtom - Genat4 10.1.p2
  double xSection = 0.;

  const double dT0 = keV;
  const double a = 20.0, b = 230.0, c = 440.0;

  /*
  double p1Z = Z*(d1 + e1*Z + f1*Z*Z);
  double p2Z = Z*(d2 + e2*Z + f2*Z*Z);
  double p3Z = Z*(d3 + e3*Z + f3*Z*Z);
  double p4Z = Z*(d4 + e4*Z + f4*Z*Z);
  */

  double p1Z = Z * (2.7965e-1 + 1.9756e-5 * Z + -3.9178e-7 * Z * Z) * barn;
  double p2Z = Z * (-1.8300e-1 + -1.0205e-2 * Z + 6.8241e-5 * Z * Z) * barn;
  double p3Z = Z * (6.7527 + -7.3913e-2 * Z + 6.0480e-5 * Z * Z) * barn;
  double p4Z = Z * (-1.9798e+1 + 2.7079e-2 * Z + 3.0274e-4 * Z * Z) * barn;

  double T0 = 15.0 * keV;
  if (Z < 1.5) {
    T0 = 40.0 * keV;
  }

  double X = math::Max(gammaEnergy, T0) / electron_mass_c2;
  xSection =
      p1Z * math::Log(1. + 2. * X) / X + (p2Z + p3Z * X + p4Z * X * X) / (1. + a * X + b * X * X + c * X * X * X);

  //  modification for low energy. (special case for Hydrogen)
  if (gammaEnergy < T0) {
    X = (T0 + dT0) / electron_mass_c2;
    double sigma =
        p1Z * math::Log(1. + 2 * X) / X + (p2Z + p3Z * X + p4Z * X * X) / (1. + a * X + b * X * X + c * X * X * X);
    double c1 = -T0 * (sigma - xSection) / (xSection * dT0);
    double c2 = 0.150;
    if (Z > 1.5) {
      c2 = 0.375 - 0.0556 * math::Log(1.0 * Z);
    }
    double y = math::Log(gammaEnergy / T0);
    xSection *= math::Exp(-y * (c1 + c2 * y));
  }
  if (xSection < 0.0) {
    xSection = 0.0;
  }
  return xSection;
}

VECCORE_ATT_HOST_DEVICE void ComptonKleinNishina::SampleByCompositionRejection(int /*Z*/, // not used
                                                                                double energyIn, double &energyOut,
                                                                                double &sinTheta)
{
  double epsilon, epsilonsq, onecost, sint2, greject;

  double E0_m = energyIn * inv_electron_mass_c2;
  double eps0 = 1. / (1. + 2. * E0_m);
  double epsilon0sq = eps0 * eps0;
  double alpha1 = -math::Log(eps0);
  double alpha2 = 0.5 * (1. - epsilon0sq);

  do {
    if (alpha1 / (alpha1 + alpha2) > UniformRandom<double>(fRandomState, fThreadId)) {
      epsilon = math::Exp(-alpha1 * UniformRandom<double>(fRandomState, fThreadId));
      epsilonsq = epsilon * epsilon;
    }
    else {
      epsilonsq = epsilon0sq + (1. - epsilon0sq) * UniformRandom<double>(fRandomState, fThreadId);
      epsilon = math::Sqrt(epsilonsq);
    }

    onecost = (1. - epsilon) / (epsilon * E0_m);
    sint2 = onecost * (2. - onecost);
    greject = 1. - epsilon * sint2 / (1. + epsilonsq);

  } while (greject < UniformRandom<double>(fRandomState, fThreadId));

  energyOut = epsilon * energyIn;
  sinTheta = (sint2 < 0.0) ? 0.0 : math::Sqrt(sint2);
}

} // end namespace impl
} // end namespace vecphys
