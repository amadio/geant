#ifndef PARTICLEGUN_H
#define PARTICLEGUN_H

#include "PrimaryGenerator.h"
#include "base/RNG.h"

#include "Geant/Typedefs.h"

using PrimaryGenerator = Geant::PrimaryGenerator;
using GeantEventInfo   = Geant::GeantEventInfo;

namespace geantphysics {

/**
 * @brief   A simple particle gun for some GeantV test applications.
 * @class   ParticleGun
 * @author  M Novak, A Ribon
 * @date    april 2016
 */


class ParticleGun : public PrimaryGenerator {
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
//  Particle_t *fPartPDG; // fca particles

  double fMass;   // rest mass of the primary [GeV]
  double fCharge; // charge of the primary
  double fPTotal; // total momentum of the primary [GeV]
  double fETotal; // total energy of the primary [GeV]

  //int fNumberoftracks; // Number of generated tracks

//#ifdef USE_ROOT
//  TRandom *fRndgen; // Random number generator
//#else
  vecgeom::RNG *fRndgen; // Random number generator from vecgeom
//#endif

public:
  ParticleGun();
  ParticleGun(int aver, int gvcode, double partekin, double xpos, double ypos, double zpos, double xdir, double ydir,
               double zdir);

  ~ParticleGun();

  // set one GeantTrack primary track properties
  virtual void InitPrimaryGenerator();
  virtual GeantEventInfo NextEvent();
  virtual void GetTrack(int n, Geant::GeantTrack &gtrack);

private:
  ParticleGun(const ParticleGun &);            // no imp.
  ParticleGun &operator=(const ParticleGun &); // no imp.

};

}

#endif
