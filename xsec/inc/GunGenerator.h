#ifndef GunGenerator_h
#define GunGenerator_h

#include "TPartIndex.h"
#include "PrimaryGenerator.h"
#ifdef USE_ROOT
#include "TRandom.h"
#else
#include "base/RNG.h"
using vecgeom::RNG;
#endif

#include "Geant/Typedefs.h"

class GunGenerator : public PrimaryGenerator {
private:
  int average;

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

  int numberoftracks;

#ifdef USE_ROOT
  TRandom *rndgen;
#else
  RNG *rndgen;
#endif

public:
  GunGenerator();
  GunGenerator(int aver, int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir, double ydir,
               double zdir);

  ~GunGenerator();

  // set one GeantTrack primary track properties
  virtual void InitPrimaryGenerator();
  virtual int NextEvent();
  virtual void GetTrack(int n, Geant::GeantTrack &gtrack);

private:
  GunGenerator(const GunGenerator &);            // no imp.
  GunGenerator &operator=(const GunGenerator &); // no imp.

  ClassDef(GunGenerator, 1)
};

#endif
