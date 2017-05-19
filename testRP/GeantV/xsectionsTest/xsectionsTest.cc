#include "GlauberGribovTotalXsc.h"
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


int main(int argc, char * argv[])
{

  Proton* p = geantphysics::Proton::Definition();
  Neutron* n = geantphysics::Neutron::Definition();
     
  geantphysics::GlauberGribovTotalXsc hnxs;

  std::ofstream writef("sumXsc.dat", std::ios::out ) ;


  writef.setf( std::ios::scientific, std::ios::floatfield );

  int i, iMax = 105;
  double kinEnergy = 10.* geant::MeV;
  int particlePDG = 2112;
  int Z = 82;
  int N = 125;
  double xs = 0;
  
  writef << iMax << std::endl;
   
  for(i=0;i<iMax;i++)
    {
      xs = hnxs.GetIsotopeCrossSection(geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetInternalCode(),
				       kinEnergy, geantphysics::Particle::GetParticleByPDGCode(particlePDG)->GetPDGMass(), Z, N) / geant::millibarn;
      
      std::cout <<"energyKin " << kinEnergy <<  " total " << xs<< " elastic " << " inelastic " <<  std::endl;
      
      writef << kinEnergy/geant::MeV <<"  \t" << 0 << "  \t"
	     << 0 << "  \t"<< 0 << "  \t"<< xs << std::endl;     
      
      kinEnergy *= 1.138;
    }
  
  exit(0);
}
