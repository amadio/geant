#include <iostream>

#include "base/VecPhys.h"

#include "GUAliasSampler.h"
#include "GUAliasTable.h"
#include "PhotoElectronSauterGavrila.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
PhotoElectronSauterGavrila::PhotoElectronSauterGavrila(Random_t *states, int tid)
    : EmModelBase<PhotoElectronSauterGavrila>(states, tid)
{
  //  Initialization();
}

VECCORE_ATT_HOST_DEVICE
PhotoElectronSauterGavrila::PhotoElectronSauterGavrila(Random_t *states, int tid, GUAliasSampler *sampler)
    : EmModelBase<PhotoElectronSauterGavrila>(states, tid, sampler)
{
}

VECCORE_ATT_HOST void PhotoElectronSauterGavrila::Initialization()
{
  if (fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, fLowEnergyLimit, fHighEnergyLimit, 100, 200);
    // note: if (Egamma/electron_mass_c2 > 50), math::Cos(theta) = 1 in Geant4

    BuildAliasTable();
  }
}

VECCORE_ATT_HOST void PhotoElectronSauterGavrila::BuildCrossSectionTablePerAtom(int /*Z*/)
{
  ; // dummy for now
}

VECCORE_ATT_HOST void PhotoElectronSauterGavrila::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (KleinNishina pdf) in the
  // input energy randge [fAliasSampler->GetIncomingMin(),fAliasSampler->GetIncomingMax()]
  //       with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number)
  // output :  p    (probability distribution) with the array dimension
  //                [fAliasSampler->GetNumEntries()][fAliasSampler->GetSamplesPerEntry()]
  //
  // storing/retrieving convention for irow and icol : p[irow x ncol + icol]

  const int nrow = fAliasSampler->GetNumEntries();
  const int ncol = fAliasSampler->GetSamplesPerEntry();

  double logxmin = math::Log(fAliasSampler->GetIncomingMin());
  double dx = (math::Log(fAliasSampler->GetIncomingMax()) - logxmin) / nrow;

  const int nintegral = 10; // temporary

  for (int i = 0; i <= nrow; ++i) {
    // for each input energy bin
    double x = math::Exp(logxmin + dx * i);

    const double ymin = -1.0;
    const double dy = 2. / ncol;
    //    const double yo = ymin + 0.5*dy;

    double sum = 0.;

    for (int j = 0; j < ncol; ++j) {

      // for each input math::Cos(theta) bin
      double normal = 0;

      // cross section weighted bin position
      for (int k = 0; k < nintegral; ++k) {

        double y = ymin + dy * (j + (0.5 + k) / nintegral);
        double fxsec = CalculateDiffCrossSectionK(Z, x, y);
        normal += fxsec;
      }
      double xsec = normal / nintegral;

      // for each output energy bin
      //      double y = yo + dy*j;
      //      double xsec = CalculateDiffCrossSectionK(Z,x,y);

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

// function implementing the angular distribution of photoelectrons

VECCORE_ATT_HOST_DEVICE double PhotoElectronSauterGavrila::CalculateDiffCrossSectionK(int /*Zelement*/, double energy,
                                                                                       double cosTheta) const
{
  // based on Geant4 : G4SauterGavrilaAngularDistribution
  // input  : energy   (incomming photon energy)
  //          cosTheta (cons(theta) of photo-electron)
  // output : dsigmaK  (differential cross section, K-shell only)

  double tau = energy / electron_mass_c2;

  double g = tau + 1.0;
  double invgamma = 1.0 / (tau + 1.0);
  double beta = math::Sqrt(tau * (tau + 2.0)) * invgamma;

  double z = 1 - beta * cosTheta;
  double z2 = z * z;
  double z4 = z2 * z2;
  double y = 1 - cosTheta * cosTheta;

  double dsigmaK = (y / z4) * (1 + 0.5 * g * (g - 1) * (g - 2) * z);

  return dsigmaK;
}

VECCORE_ATT_HOST_DEVICE double PhotoElectronSauterGavrila::CalculateDiffCrossSection(int Zelement, double energy,
                                                                                      double cosTheta) const
{
  // based on Geant4 : G4SauterGavrilaAngularDistribution
  // input  : energy   (incomming photon energy)
  //          cosTheta (cons(theta) of photo-electron)
  // output : dsigma  (differential cross section, K-shell + L-shells)

  double tau = energy / electron_mass_c2;

  double g = tau + 1.0;
  double invgamma = 1.0 / (tau + 1.0);
  double beta = math::Sqrt(tau * (tau + 2.0)) * invgamma;

  double g2 = g * g;
  double g3 = g2 * g;
  double g4 = g2 * g2;
  //  double g5 = g2*g3;

  double term = math::Log(g * (1. + beta)) / (g * beta);

  double sigmaL2 = g3 - 5. * g2 + 24. * g - 16. + (g2 + 3 * g - 8) * term;
  double sigmaL3 = 4. * g3 - 6. * g2 + 5. * g + 3. + (g2 - 3 * g + 4) * term;

  double sigma23 = sigmaL2 + sigmaL3;

  double JK = 125. / Zelement + 3.5;
  double JL = 1.2;

  double R2 = sigmaL2 / sigma23;
  double R3 = sigmaL3 / sigma23;

  double PK = 1 - 1. / JK;
  double PL = (1. - PK) * (1 - 1. / JL);
  double PL2 = (1. - PK - PL) * R2;
  double PL3 = (1. - PK - PL) * R3;

  double z = 1 - beta * cosTheta;
  double z2 = z * z;
  double z3 = z2 * z;
  double z4 = z2 * z2;
  double z5 = z2 * z3;
  double y = 1 - cosTheta * cosTheta;

  double dsigmaK = (y / z4) * (1 + 0.5 * g * (g - 1) * (g - 2) * z) * PK;
  double dsigmaL1 = dsigmaK;

  double coeff = math::Sqrt((g + 1) * tau) / math::Pow(g * tau, 5.0);

  double dsigmaL2 = g * (3. * g + 1) / (2 * z4) - g2 * (9 * g2 + 30 * g - 7) / (8 * z3) +
                    g3 * (g3 + 6 * g2 + 11 * g - 2) / (4 * z2) - g4 * (tau * (g + 7)) / (8 * z) +
                    y * ((g + 1) / z5 - g * (g + 1) / z4 + g2 * (3 * g + 1) * (g2 - 1) / z3);

  dsigmaL2 *= coeff;

  double dsigmaL3 = g * (1 - 3 * g) / (2 * z4) + g2 * (3 * g2 - 1) / z3 + g3 * (g3 - 3 * g2 + 2 * g + 1) / z3 -
                    g4 * (g - 2) * (g - 1) / (2 * z) +
                    y * ((g + 1) / z5 - g * (g + 1) * (2 * g - 1) / z4 + g2 * (g2 - 1) / z3);

  dsigmaL3 *= coeff;

  double dsigma = PK * dsigmaK + PL * dsigmaL1 + PL2 * dsigmaL2 + PL3 * dsigmaL3;
  return dsigma;
}

VECCORE_ATT_HOST double PhotoElectronSauterGavrila::GetG4CrossSection(int Z, double energy)
{
  // G4PEEffectFluoModel::ComputeCrossSectionPerAtom

  // This method may be used only if G4MaterialCutsCouple pointer
  //   has been set properly

  // Sandia parameterization for Z < 100
  //  int Z = zElement;

  int fCumulInterval[101] = {0};
  double fSandiaCof[4] = {0.0};

  fCumulInterval[0] = 1;

  // scan - move it to constructor or use pre-built table
  for (int iz = 1; iz < 101; ++iz) {
    fCumulInterval[iz] = fCumulInterval[iz - 1] + fNbOfIntervals[iz];
  }

  double Emin = fSandiaTable[fCumulInterval[Z - 1]][0] * keV;

  int interval = fNbOfIntervals[Z] - 1;
  int row = fCumulInterval[Z - 1] + interval;

  while ((interval > 0) && (energy < fSandiaTable[row][0] * keV)) {
    --interval;
    row = fCumulInterval[Z - 1] + interval;
  }

  if (energy >= Emin) {
    double AoverAvo = Z * amu / fZtoAratio[Z];
    fSandiaCof[0] = AoverAvo * funitc[1] * fSandiaTable[row][1];
    fSandiaCof[1] = AoverAvo * funitc[2] * fSandiaTable[row][2];
    fSandiaCof[2] = AoverAvo * funitc[3] * fSandiaTable[row][3];
    fSandiaCof[3] = AoverAvo * funitc[4] * fSandiaTable[row][4];
  }
  else {
    fSandiaCof[0] = fSandiaCof[1] = fSandiaCof[2] = fSandiaCof[3] = 0.;
  }

  //   CurrentCouple()->GetMaterial()
  //     ->GetSandiaTable()->GetSandiaCofPerAtom((int)Z, energy, fSandiaCof);

  double energy2 = energy * energy;
  double energy3 = energy * energy2;
  double energy4 = energy2 * energy2;

  return fSandiaCof[0] / energy + fSandiaCof[1] / energy2 + fSandiaCof[2] / energy3 + fSandiaCof[3] / energy4;

  double xSection = 1.0 * 10e-24;
  // dummy for now
  return xSection;
}

VECCORE_ATT_HOST_DEVICE void PhotoElectronSauterGavrila::SampleByCompositionRejection(int Z, // not used
                                                                                       double energyIn,
                                                                                       double &energyOut, double &sint)
{
  // use the scalar implementation which is equivalent to Geant4
  energyOut = GetPhotoElectronEnergy<ScalarBackend>(energyIn, Z);

  // sample angular direction according to SauterGavrilaAngularDistribution

  double tau = energyIn / electron_mass_c2;
  //  static
  const double taulimit = 50.0;

  double cost = -1.0;

  if (tau > taulimit) {
    cost = 1.0; // set to the primary direction
  }
  else {
    // algorithm according Penelope 2008 manual and
    // F.Sauter Ann. Physik 9, 217(1931); 11, 454(1931).

    double gamma = tau + 1.;
    double beta = math::Sqrt(tau * (tau + 2.)) / gamma;
    double A = (1 - beta) / beta;
    double Ap2 = A + 2.;
    double B = 0.5 * beta * gamma * (gamma - 1.) * (gamma - 2.);
    double grej = 2. * (1. + A * B) / A;
    double z, g;
    do {
      double q = UniformRandom<double>(fRandomState, fThreadId);
      z = 2 * A * (2 * q + Ap2 * math::Sqrt(q)) / (Ap2 * Ap2 - 4 * q);
      g = (2 - z) * (1.0 / (A + z) + B);

    } while (g < UniformRandom<double>(fRandomState, fThreadId) * grej);

    cost = 1 - z;
  }

  sint = math::Sqrt((1 + cost) * (1 - cost));
}

} // end namespace impl
} // end namespace vecphys
