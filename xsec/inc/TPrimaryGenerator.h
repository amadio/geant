#ifndef TPrimaryGenerator_H
#define TPrimaryGenerator_H

#include "TPartIndex.h"
#include "Geant/Typedefs.h"
#include "GeantFwd.h"
#ifndef GEANT_NVCC
#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#endif

class TPrimaryGenerator {
private:
  int fPDG;         // PDG code of parimary particles
  double fPartEkin; // kinetic energy of the primary [GeV]
  double fXPos;     // (x,y,z) position of the primary particles
  double fYPos;
  double fZPos;
  double fXDir; // direction vector of the primary particles
  double fYDir;
  double fZDir;
  // additional members
  int fGVPartIndex; // GV particle index of the primary
  Particle_t *fPartPDG;
  double fMass;   // rest mass of the primary [GeV]
  double fCharge; // charge of the primary
  double fPTotal; // total momentum of the primary [GeV]
  double fETotal; // total energy of the primary [GeV]

public:
  TPrimaryGenerator();
  TPrimaryGenerator(int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir, double ydir,
                    double zdir);
  virtual ~TPrimaryGenerator();

  void SetEkin(double partekin) { fPartEkin = partekin; }
  double GetEkin() const { return fPartEkin; }
  void SetParticleByPDGCode(int pdgcode);
  int GetParticlePDGCode() const { return fPDG; }
  void SetParticleXYZPosition(double x, double y, double z) {
    fXPos = x;
    fYPos = y;
    fZPos = z;
  }
  double GetParticleXPos() const { return fXPos; }
  double GetParticleYPos() const { return fYPos; }
  double GetParticleZPos() const { return fZPos; }
  void SetParticleXYZDir(double xdir, double ydir, double zdir);
  double GetParticleXDir() const { return fXDir; }
  double GetParticleYDir() const { return fYDir; }
  double GetParticleZDir() const { return fZDir; }

  // --
  int GetParticleGVIndex() const { return fGVPartIndex; }
  const Particle_t *GetParticlePDG() const { return fPartPDG; }
  double GetParticleMass() const { return fMass; }
  double GetParticleCharge() const { return fCharge; }
  double GetparticlePTotal() const { return fPTotal; }
  double GetparticleETotal() const { return fETotal; }

  // set one GeantTrack primary track properties
  void InitPrimaryTrack(Geant::GeantTrack &gtrack);

private:
  void InitPrimaryGenerator();

  TPrimaryGenerator(const TPrimaryGenerator &);            // no imp.
  TPrimaryGenerator &operator=(const TPrimaryGenerator &); // no imp.

#ifdef USE_ROOT
  ClassDefNV(TPrimaryGenerator, 1)
#endif
};

#endif
