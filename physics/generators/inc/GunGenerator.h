#ifndef GunGenerator_h
#define GunGenerator_h

#include "TPartIndex.h"
#include "Geant/PrimaryGenerator.h"
#ifdef USE_ROOT
#include "TRandom.h"
#else
#include "base/RNG.h"
using vecgeom::RNG;
#endif

#include "Geant/Typedefs.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class TaskData;

class GunGenerator : public PrimaryGenerator {
private:
  int fAverage; // Average number of tracks for Poisson distribution

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

#ifdef USE_ROOT
  TRandom *fRndgen; // Random number generator
#else
  RNG *fRndgen; // Random number generator
#endif

public:
  GunGenerator();
  GunGenerator(int aver, int partpdg, double partekin, double xpos, double ypos, double zpos, double xdir, double ydir,
               double zdir);

  ~GunGenerator();

  // set one Track primary track properties
  virtual void InitPrimaryGenerator();
  virtual EventInfo NextEvent(geant::TaskData* td);
  virtual void GetTrack(int n, geant::Track &gtrack, geant::TaskData* td);

private:
  GunGenerator(const GunGenerator &);            // no imp.
  GunGenerator &operator=(const GunGenerator &); // no imp.
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
