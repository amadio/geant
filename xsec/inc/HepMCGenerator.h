
#ifndef HepMCGenerator_h
#define HepMCGenerator_h

#include "TPartIndex.h"
#include "PrimaryGenerator.h"

#if __cplusplus >= 201103L
#include "HepMC/Reader.h"
#include "HepMC/GenEvent.h"
#include "HepMC/Search/FindParticles.h"
#endif

#ifdef USE_VECGEOM_NAVIGATOR
#include "volumes/Particle.h"
using vecgeom::Particle;
#else
class TParticlePDG;
#endif

class HepMCGenerator : public PrimaryGenerator {
private:
#if __cplusplus >= 201103L
  HepMC::Reader *input_file;
  HepMC::FindParticles *search;
#endif

public:
  HepMCGenerator();
  HepMCGenerator(std::string &filename);
  ~HepMCGenerator();

  // set one GeantTrack primary track properties
  virtual void InitPrimaryGenerator();
  virtual Int_t NextEvent();
  virtual void GetTrack(Int_t n, Geant::GeantTrack &gtrack);
  // used from Geant4 test-complex to take one primary track
  void GetTrack(Int_t n, Double_t &tpx, Double_t &tpy, Double_t &tpz, Double_t &te, Double_t &x0, Double_t &y0,
                Double_t &z0, Int_t &pdg);

private:
  HepMCGenerator(const HepMCGenerator &); // no imp.
  HepMCGenerator &operator=(const HepMCGenerator &); // no imp.

  ClassDef(HepMCGenerator, 1)
};

#endif
