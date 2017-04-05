
#include "Particle.h"

//#include "PhysicalConstants.h"

namespace geantphysics {

std::vector<Particle*> Particle::gTheParticleTable;
std::vector<Particle*> Particle::gInternalParticleCodes;

Particle::Particle(const std::string &name, int pdgcode, int intcode, double charge, double mass)
: fName(name), fIndex(-1), fInternalCode(intcode), fPDGCode(pdgcode), fPDGCharge(charge), fPDGMass(mass) {
  fIndex = gTheParticleTable.size();
  gTheParticleTable.push_back(this);
  //
  unsigned long icode = intcode;
  if (gInternalParticleCodes.size()<icode+1) {
    gInternalParticleCodes.resize(icode+1,nullptr);
  }
  gInternalParticleCodes[icode] = this;
}

} // namespace geantphysics
