
#include "ELossTableRegister.h"

namespace geantphysics {

ELossTableRegister& ELossTableRegister::Instance() {
  static ELossTableRegister instance;
  return instance;
}


void ELossTableRegister::RegisterEnergyLossProcessForParticle(int internalpartindx, EMPhysicsProcess *elossprocess) {
  unsigned long ipartindx = internalpartindx;
  if (fELossProcessesPerParticle.size()<ipartindx+1) {
    fELossProcessesPerParticle.resize(ipartindx+1);
  }
  fELossProcessesPerParticle[ipartindx].push_back(elossprocess);
}


void ELossTableRegister::Clear() {
  for (unsigned long i=0; i<fELossProcessesPerParticle.size(); ++i) {
    fELossProcessesPerParticle[i].clear();
  }
  fELossProcessesPerParticle.clear();
}

} // namespace geantphysics
