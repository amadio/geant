#include "GlauberGribovTotalXsc.h"
#include "GlauberGribovInelasticXsc.h"
#include "GlauberGribovElasticXsc.h"
#include "Particle.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "Proton.h"
#include "Neutron.h"
#include "Electron.h"

#include "SystemOfUnits.h"
#include "PhysicalConstants.h"

using geantphysics::Particle;
using geantphysics::Proton;
using geantphysics::Neutron;


int main(int /*argc*/, char** /*argv*/) {
  // Proton*  p = geantphysics::Proton::Definition();
  // Neutron* n = geantphysics::Neutron::Definition();

  geantphysics::GlauberGribovTotalXsc totxs;
  geantphysics::GlauberGribovInelasticXsc inexs;
  geantphysics::GlauberGribovElasticXsc elaxs;

  std::ofstream writef("sumXsc.dat", std::ios::out ) ;

  writef.setf( std::ios::scientific, std::ios::floatfield );

  int i, iMax = 105;
  double kinEnergy = 10.* geant::MeV;
  int particlePDG = 2112;
  int Z = 82;
  int N = 125;
  double txs, ixs, exs;

  writef << iMax << std::endl;
  for (i=0; i<iMax; i++) {
    txs = totxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    ixs = inexs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    exs = elaxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
			       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;

    std::cout <<"energyKin " << kinEnergy <<  " total " << txs << " elastic " << exs << " inelastic " << ixs <<  std::endl;

    writef << kinEnergy/geant::MeV <<"  \t" << 0 << "  \t"
     << 0 << "  \t"<< 0 << "  \t"<< txs << std::endl;

    kinEnergy *= 1.138;
  }

  return 0;
}
