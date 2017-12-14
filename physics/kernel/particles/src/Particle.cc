
#include "Particle.h"

//#include "PhysicalConstants.h"

namespace geantphysics {

std::vector<Particle*> Particle::gTheParticleTable;
std::vector<Particle*> Particle::gInternalParticleCodes;
std::map<unsigned int, unsigned int>   Particle::gPDGtoInternalCode;
std::map<const std::string, Particle*> Particle::gNametoParticle;

Particle::Particle(const std::string &name, int pdgcode, int intcode, double mass, double charge)
: fName(name), fIndex(-1), fInternalCode(intcode), fPDGCode(pdgcode), fPDGMass(mass), fPDGCharge(charge) {
  fIndex = gTheParticleTable.size();
  gTheParticleTable.push_back(this);
  //
  unsigned long icode = intcode;
  if (gInternalParticleCodes.size()<icode+1) {
    gInternalParticleCodes.resize(icode+1,nullptr);
  }
  gInternalParticleCodes[icode] = this;
  //
  gPDGtoInternalCode[pdgcode] = icode;
  gNametoParticle[name]       = this;
}

} // namespace geantphysics
