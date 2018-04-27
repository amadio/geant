
#include "Geant/CutConverterForPositron.h"

#include "Geant/Material.h"
#include "Geant/MaterialProperties.h"
#include "Geant/Element.h"

#include <cmath>
#include <iostream>
#include "Geant/math_wrappers.h"

namespace geantphysics {

CutConverterForPositron::CutConverterForPositron(int numebins, double mincutenergy, double maxcutenergy)
    : CutConverter(2, numebins, mincutenergy, maxcutenergy)
{
  if (fMinCutEnergy >= fMaxCutEnergy) {
    std::cerr << "  *** ERROR in CutConverterForPositron::CutConverterForPositron() " << std::endl
              << "       minenergy = " << mincutenergy / geant::units::GeV
              << " [GeV] >= maxenergy = " << maxcutenergy / geant::units::GeV << " [GeV]" << std::endl;
    exit(-1);
  }
  Initialise();
}

// must be called before using the Convert method if new element has been inserted into the Element table!
void CutConverterForPositron::Initialise()
{
  CutConverter::Initialise();
  BuildElossOrAbsXsecTable();
}

CutConverterForPositron::~CutConverterForPositron()
{
}

double CutConverterForPositron::ComputeELossOrAbsXsecPerAtom(double zet, double ekin)
{
  const double cbr1  = 0.02;
  const double cbr2  = -5.7e-5;
  const double cbr3  = 1.;
  const double cbr4  = 0.072;
  const double tlow  = 10. * geant::units::keV;
  const double thigh = 1. * geant::units::GeV;
  const double mass  = geant::units::kElectronMassC2;
  const double taul  = tlow / mass;
  const double cpot  = 1.6e-5 * geant::units::MeV;
  const double fact  = geant::units::kTwoPi * geant::units::kElectronMassC2 * geant::units::kClassicElectronRadius *
                      geant::units::kClassicElectronRadius;

  double ionpot    = cpot * Math::Exp(0.9 * Math::Log(zet)) / mass;
  double ionpotlog = Math::Log(ionpot);

  // calculate approximated dE/dx for electrons
  double tau  = ekin / mass;
  double dEdx = 0.0;
  if (tau < taul) {
    double t1    = taul + 1.;
    double t2    = taul + 2.;
    double tsq   = taul * taul;
    double beta2 = taul * t2 / (t1 * t1);
    // this is different compared to e-
    double f = 2. * Math::Log(taul) -
               (6. * taul + 1.5 * tsq - taul * (1. - tsq / 3.) / t2 - tsq * (0.5 - tsq / 12.) / (t2 * t2)) / (t1 * t1);
    dEdx        = (Math::Log(2. * taul + 4.) - 2. * ionpotlog + f) / beta2;
    dEdx        = fact * zet * dEdx;
    double clow = dEdx * std::sqrt(taul);
    dEdx        = clow / std::sqrt(tau);
  } else {
    double t1    = tau + 1.;
    double t2    = tau + 2.;
    double tsq   = tau * tau;
    double beta2 = tau * t2 / (t1 * t1);
    // this is different compared to e-
    double f = 2. * Math::Log(tau) -
               (6. * tau + 1.5 * tsq - tau * (1. - tsq / 3.) / t2 - tsq * (0.5 - tsq / 12.) / (t2 * t2)) / (t1 * t1);
    dEdx = (Math::Log(2. * tau + 4.) - 2. * ionpotlog + f) / beta2;
    dEdx = fact * zet * dEdx;
    // loss from bremsstrahlung follows
    double cbrem = (cbr1 + cbr2 * zet) * (cbr3 + cbr4 * Math::Log(ekin / thigh));
    cbrem        = 0.1 * zet * (zet + 1.) * cbrem * tau / beta2;
    dEdx += fact * cbrem;
  }
  return dEdx;
}

} // namespace geantphysics
