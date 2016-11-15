#include "IonisationMoller.h"
#include "GUAliasSampler.h"
#include "GUAliasTable.h"

#include "base/VecPhys.h"

namespace vecphys {
inline namespace VECPHYS_IMPL_NAMESPACE {

VECCORE_ATT_HOST
IonisationMoller::IonisationMoller(Random_t *states, int tid) : EmModelBase<IonisationMoller>(states, tid)
{
  fDeltaRayThreshold = 1.0 * keV; // temporary: should be set from a cut table
  SetLowEnergyLimit(0.1 * keV);

  Initialization();
}

VECCORE_ATT_HOST_DEVICE
IonisationMoller::IonisationMoller(Random_t *states, int tid, GUAliasSampler *sampler)
    : EmModelBase<IonisationMoller>(states, tid, sampler)
{
  fDeltaRayThreshold = 1.0 * keV; // temporary: should be set from a cut table
  SetLowEnergyLimit(0.1 * keV);
}

VECCORE_ATT_HOST void IonisationMoller::Initialization()
{
  if (fSampleType == kAlias) {
    fAliasSampler = new GUAliasSampler(fRandomState, fThreadId, 1.e-4, 1.e+6, 100, 100);
    BuildAliasTable();
  }
}

VECCORE_ATT_HOST void IonisationMoller::BuildCrossSectionTablePerAtom(int /*Z*/)
{
  ; // dummy for now
}

VECCORE_ATT_HOST void IonisationMoller::BuildPdfTable(int Z, double *p)
{
  // Build the probability density function (MollerBhabha pdf) in the
  // input energy randge [fMinX,fMaxX] with an equal logarithmic bin size
  //
  // input  :  Z    (atomic number) - not used for the atomic independent model
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

    // e-e- (Moller) only for now
    double ymin = fDeltaRayThreshold;      // minimum delta-ray energy
    double dy   = (x / 2.0 - ymin) / ncol; // maximum x/2.0
    double yo   = ymin + 0.5 * dy;

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

// function implementing the cross section for MollerBhabha

VECCORE_ATT_HOST_DEVICE double IonisationMoller::CalculateDiffCrossSection(int /*Zelement*/, double kineticEnergy,
                                                                           double deltaRayEnergy) const
{
  // based on Geant3 : Simulation of the delta-ray production (PHY331-1)
  // input  : kineticEnergy (incomming photon energy)
  //          deltaRayEnergy (scattered photon energy)
  // output : dcross (differential cross section)

  double dcross = 0.0;

  double tau    = kineticEnergy / electron_mass_c2;
  double gam    = tau + 1.0;
  double gamma2 = gam * gam;

  double epsil = deltaRayEnergy / kineticEnergy;

  // Moller (e-e-) scattering only
  // Bhabha (e+e-) scattering not implemented for now

  double fgam = (2 * gam - 1.0) / gamma2;
  double x    = 1 / epsil;
  double y    = 1 / (1.0 - epsil);

  dcross = 1 - fgam + x * (x - fgam) + y * (y - fgam);

  return dcross;
}

VECCORE_ATT_HOST double IonisationMoller::GetG4CrossSection(int Z, double kineticEnergy)
{
  double cross = 0.0;

  // temporary - set by material
  double cutEnergy = fDeltaRayThreshold;
  double maxEnergy = 1.0 * TeV;

  double tmax = 0.5 * kineticEnergy;
  tmax        = math::Min(maxEnergy, tmax);

  if (cutEnergy < tmax) {
    double xmin   = cutEnergy / kineticEnergy;
    double xmax   = tmax / kineticEnergy;
    double tau    = kineticEnergy / electron_mass_c2;
    double gam    = tau + 1.0;
    double gamma2 = gam * gam;
    double beta2  = tau * (tau + 2) / gamma2;

    // Moller (e-e-) scattering

    double gg = (2.0 * gam - 1.0) / gamma2;
    cross     = ((xmax - xmin) * (1.0 - gg + 1.0 / (xmin * xmax) + 1.0 / ((1.0 - xmin) * (1.0 - xmax))) -
             gg * math::Log(xmax * (1.0 - xmin) / (xmin * (1.0 - xmax)))) /
            beta2;
  }
  cross *= Z * twopi_mc2_rcl2 / kineticEnergy;

  return cross;
}

VECCORE_ATT_HOST_DEVICE void IonisationMoller::SampleByCompositionRejection(int /*Z*/, // not used
                                                                            double kineticEnergy,
                                                                            double &deltaKinEnergy, double &sinTheta)
{
  // temporary - set by material
  double cutEnergy = fDeltaRayThreshold;
  double maxEnergy = 1.0 * TeV;

  // based on G4MollerBhabhaModel::SampleSecondaries

  double tmin = cutEnergy;
  double tmax = 0.5 * kineticEnergy;

  if (maxEnergy < tmax) {
    tmax = maxEnergy;
  }

  if (tmin >= tmax) {
    return;
  }

  double energy = kineticEnergy + electron_mass_c2;
  double xmin   = tmin / kineticEnergy;
  double xmax   = tmax / kineticEnergy;
  double gam    = energy / electron_mass_c2;
  double gamma2 = gam * gam;
  //  double beta2  = 1.0 - 1.0/gamma2;
  double x, z, q, grej;

  // Moller (e-e-) scattering
  double gg = (2.0 * gam - 1.0) / gamma2;
  double y  = 1.0 - xmax;
  grej      = 1.0 - gg * xmax + xmax * xmax * (1.0 - gg + (1.0 - gg * y) / (y * y));

  do {
    q = UniformRandom<double>(fRandomState, fThreadId);
    x = xmin * xmax / (xmin * (1.0 - q) + xmax * q);
    y = 1.0 - x;
    z = 1.0 - gg * x + x * x * (1.0 - gg + (1.0 - gg * y) / (y * y));
  } while (grej * UniformRandom<double>(fRandomState, fThreadId) > z);

  deltaKinEnergy = x * kineticEnergy;

  double totalMomentum = math::Sqrt(kineticEnergy * (kineticEnergy + 2.0 * electron_mass_c2));

  double deltaMomentum = math::Sqrt(deltaKinEnergy * (deltaKinEnergy + 2.0 * electron_mass_c2));
  double cost          = deltaKinEnergy * (energy + electron_mass_c2) / (deltaMomentum * totalMomentum);
  if (cost > 1.0) {
    cost = 1.0;
  }
  double sint2 = (1.0 - cost) * (1.0 + cost);

  sinTheta = (sint2 < 0.0) ? 0.0 : math::Sqrt(sint2);
}

} // end namespace impl
} // end namespace vecphys
