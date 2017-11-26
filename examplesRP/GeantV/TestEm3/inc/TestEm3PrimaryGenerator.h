
#ifndef TESTEM3PRIMARYGENERATOR_H
#define TESTEM3PRIMARYGENERATOR_H

// from GeantV
#include "PrimaryGenerator.h"
// for GEANT_IMPL_NAMESPACE
#include "Geant/Typedefs.h"

namespace GEANT_IMPL_NAMESPACE {
  namespace Geant {
    class GeantTrack;
    class GeantTaskData;
    class GeantEventInfo;
  }
}

// geantphysics classdef
namespace geantphysics {
  class Particle;
}

#include <string>

namespace userapplication {

/**
 * @brief User primary generator for TestEm3App.
 *
 * The primary generator is a simple particle gun with configurable primary particle type and kinetic energy.
 *
 * @class   TestEm3PrimaryGenerator
 * @author  M Novak
 * @date    November 2017
 */

class TestEm3DetectorConstruction;

class TestEm3PrimaryGenerator : public Geant::PrimaryGenerator {
public:
  // CTR DTR
  TestEm3PrimaryGenerator(TestEm3DetectorConstruction *det);
 ~TestEm3PrimaryGenerator();

  // public setters/getters
  void SetPrimaryParticleName(std::string pname) { fPrimaryParticleName = pname; }
  const std::string & GetPrimaryParticleName() const { return fPrimaryParticleName; }

  void   SetPrimaryParticleEnergy(double ekin) { fPrimaryEkin = ekin; }
  double GetPrimaryParticleEnergy() const { return fPrimaryEkin; }

  void SetNumberOfPrimaryParticlePerEvent(int val) { fPrimaryPerEvent = val; }
  int  GetNumberOfPrimaryParticlePerEvent() const { return fPrimaryPerEvent; }

  const geantphysics::Particle* GetPrimaryParticle() const { return fParticle; }

  // interface methods
  virtual void InitPrimaryGenerator();
  virtual Geant::GeantEventInfo NextEvent(Geant::GeantTaskData* td);
  virtual void GetTrack(int n, Geant::GeantTrack &gtrack, Geant::GeantTaskData* td);

private:
 TestEm3PrimaryGenerator() = delete;
 TestEm3PrimaryGenerator(const TestEm3PrimaryGenerator &) = delete;
 TestEm3PrimaryGenerator &operator=(const TestEm3PrimaryGenerator &) = delete;

private:
  std::string  fPrimaryParticleName; // name of the primary particle
  int          fPrimaryPerEvent;     // number of primary particle to be generated per event
  int    fPDG;          // PDG code of parimary particles
  int    fGVPartIndex;  // internal GV particle index of the primary
  double fPrimaryEkin;  // kinetic energy of the primary in internal [energy] unit
  double fXPos;         // (x,y,z) position of the primary particles in internal [length] unit
  double fYPos;
  double fZPos;
  double fXDir;          // direction vector of the primary particles
  double fYDir;
  double fZDir;
  //
  double fMass;         // rest mass of the primary in internal [energy] unit
  double fCharge;       // charge of the primary in internal [charge] unit
  double fETotal;       // total energy of the primary in internal [energy] unit
  double fPTotal;       // total momentum of the primary in internal [energy] unit
  //
  const geantphysics::Particle    *fParticle; // the primary particle
  TestEm3DetectorConstruction        *fDetector; // the detector
};

}       // namespace userapplication

#endif  // TESTEM3PRIMARYGENERATOR_H
