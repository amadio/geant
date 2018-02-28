#include "Geant/GlauberGribovTotalXsc.h"
#include "Geant/GlauberGribovInelasticXsc.h"
#include "Geant/GlauberGribovElasticXsc.h"
#include "Geant/Particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Geant/Proton.h"
#include "Geant/Neutron.h"
#include "Geant/KaonMinus.h"
#include "Geant/Electron.h"

#include "Geant/SystemOfUnits.h"
#include "Geant/PhysicalConstants.h"

using geantphysics::Particle;
using geantphysics::Proton;
using geantphysics::Neutron;
using geantphysics::KaonMinus;

int main(int /*argc*/, char ** /*argv*/)
{
  geantphysics::Proton::Definition();
  geantphysics::Neutron::Definition();
  geantphysics::KaonMinus::Definition();

  geantphysics::GlauberGribovTotalXsc totxs;
  geantphysics::GlauberGribovInelasticXsc inexs;
  geantphysics::GlauberGribovElasticXsc elaxs;

  std::ofstream writef("sumXsc.dat", std::ios::out);

  writef.setf(std::ios::scientific, std::ios::floatfield);

  double kinEnergy = 10. * geant::units::MeV;
  double maxEnergy = 1000 * geant::units::GeV;
  int particlePDG  = -321;
  // int particlePDG = 2212;
  int Z = 82;
  // number of nucleons
  int N = 207;
  double txs, ixs, exs;

  while (kinEnergy < maxEnergy) {

    txs = totxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
                                       kinEnergy,
                                       geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) /
          geant::units::millibarn;

    ixs = inexs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
                                       kinEnergy,
                                       geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) /
          geant::units::millibarn;

    exs = elaxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
                                       kinEnergy,
                                       geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) /
          geant::units::millibarn;

    std::cout << "energyKin " << kinEnergy << " total " << txs << " elastic " << exs << " inelastic " << ixs
              << std::endl;

    writef << kinEnergy / geant::units::MeV << "  \t" << ixs << "  \t" << exs << "  \t" << 0 << "  \t" << txs
           << std::endl;

    kinEnergy *= 1.001;
  }

  return 0;
}
