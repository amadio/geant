#ifndef Particle_H
#define Particle_H

#include "base/Global.h"
#include "Geant/Config.h"

#ifdef VECCORE_CUDA
#include "base/Map.h"
#include "base/Vector.h"
#include <string.h>
#else
#include <map>
#include <vector>
#include <string>
#endif
#include <fstream>
#include <math.h>

#include <iostream>

namespace geant {

VECGEOM_DEVICE_FORWARD_DECLARE(class Material;);
VECGEOM_DEVICE_FORWARD_DECLARE(class Particle;);

inline namespace GEANT_IMPL_NAMESPACE {

#ifdef VECCORE_CUDA
class ParticleOld;
extern VECCORE_ATT_DEVICE vecgeom::map<int, ParticleOld> *fParticlesDev; // Particle list indexed by PDG code
extern vecgeom::map<int, ParticleOld> *fParticlesHost;                           // Particle list indexed by PDG code
#endif

class ParticleOld {
public:
  class Decay;
#ifdef VECCORE_CUDA
  using Map_t         = vecgeom::map<int, ParticleOld>;
  using VectorDecay_t = vecgeom::Vector<Decay>;
  using VectorInt_t   = vecgeom::Vector<int>;
#else
  using Map_t         = std::map<int, ParticleOld>;
  using VectorDecay_t = std::vector<Decay>;
  using VectorInt_t   = std::vector<int>;
#endif

  VECCORE_ATT_HOST_DEVICE
  ParticleOld();
  VECCORE_ATT_HOST_DEVICE
  ParticleOld(const char *name, int pdg, bool matter, const char *pclass, int pcode, double charge, double mass,
           double width, int isospin, int iso3, int strange, int flavor, int track, int code = -1);

  VECCORE_ATT_HOST_DEVICE
  ParticleOld(const ParticleOld &other);

  VECCORE_ATT_HOST_DEVICE
  ParticleOld &operator=(const ParticleOld &part);

  VECCORE_ATT_HOST_DEVICE
  static void CreateParticles();
  VECCORE_ATT_HOST_DEVICE
  const char *Name() const { return fName; }
  VECCORE_ATT_HOST_DEVICE
  int PDG() const { return fPDG; }
  VECCORE_ATT_HOST_DEVICE
  bool Matter() const { return fMatter; }
  VECCORE_ATT_HOST_DEVICE
  double Mass() const { return fMass; }
  const char *Class() const { return fClass; }
  VECCORE_ATT_HOST_DEVICE
  int Pcode() const { return fPcode; }
  VECCORE_ATT_HOST_DEVICE
  double Charge() const { return fCharge; }
  VECCORE_ATT_HOST_DEVICE
  double Width() const { return fWidth; }
  VECCORE_ATT_HOST_DEVICE
  int Isospin() const { return fIsospin; }
  VECCORE_ATT_HOST_DEVICE
  int Iso3() const { return fIso3; }
  VECCORE_ATT_HOST_DEVICE
  int Strange() const { return fStrange; }
  VECCORE_ATT_HOST_DEVICE
  int Flavor() const { return fFlavor; }
  VECCORE_ATT_HOST_DEVICE
  int Track() const { return fTrack; }
  VECCORE_ATT_HOST_DEVICE
  int Ndecay() const { return fNdecay; }
  VECCORE_ATT_HOST_DEVICE
  int Code() const { return fCode; }

  VECCORE_ATT_HOST_DEVICE
  void SetCode(int code) { fCode = code; }

  const VectorDecay_t &DecayList() const { return fDecayList; }
#ifndef VECCORE_CUDA
  static void ReadFile(std::string infilename, bool outfile=false);
#endif

#ifndef VECCORE_CUDA
  static const ParticleOld &GetParticle(int pdg)
  {
    if (fParticles->find(pdg) != fParticles->end()) return (*fParticles)[pdg];
    static ParticleOld p;
    std::cout << __func__ << "::pdg:" << pdg << " does not exist" << std::endl;
    return p;
  }
#else
  static const ParticleOld &GetParticle(int pdg)
  {
    if (fParticlesHost->find(pdg) != fParticlesHost->end()) return (*fParticlesHost)[pdg];
    static ParticleOld p;
    printf(" pdg %d does not exist\n", pdg);
    return p;
  }
  VECCORE_ATT_DEVICE
  static const ParticleOld &GetParticleDev(int pdg)
  {
    if (fParticlesDev->find(pdg) != fParticlesDev->end()) return (*fParticlesDev)[pdg];
    // ParticleOld p;
    printf(" pdg %d does not exist\n", pdg);
    return (*fParticlesDev)[1];
  }
#endif

#ifndef VECCORE_CUDA
  void NormDecay();

  friend std::ostream &operator<<(std::ostream &os, const ParticleOld &part);

#endif
  VECCORE_ATT_HOST_DEVICE
  void AddDecay(const Decay &decay)
  {
    fDecayList.push_back(decay);
    fNdecay = fDecayList.size();
  }
  VECCORE_ATT_HOST_DEVICE
  static const Map_t &Particles()
  {
#ifndef VECCORE_CUDA
    return *fParticles;
#else
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
    return *fParticlesHost;
#else
    return *fParticlesDev;
#endif
#endif
  }

  class Decay {
  public:
    VECCORE_ATT_HOST_DEVICE
    Decay() : fType(0), fBr(0) {}
    VECCORE_ATT_HOST_DEVICE
    Decay(const Decay &other) : fType(other.fType), fBr(other.fBr), fDaughters(other.fDaughters) {}
    VECCORE_ATT_HOST_DEVICE
    Decay(int type, double br, const VectorInt_t &daughters) : fType(type), fBr(br), fDaughters(daughters) {}
    VECCORE_ATT_HOST_DEVICE
    void Clear()
    {
      fType = 0;
      fBr   = 0;
      fDaughters.clear();
    }

    VECCORE_ATT_HOST_DEVICE
    int Type() const { return fType; }
    VECCORE_ATT_HOST_DEVICE
    double Br() const { return fBr; }
    VECCORE_ATT_HOST_DEVICE
    const VectorInt_t &Daughters() const { return fDaughters; }
    VECCORE_ATT_HOST_DEVICE
    int NDaughters() const { return fDaughters.size(); }
    VECCORE_ATT_HOST_DEVICE
    int Daughter(int i) const { return fDaughters[i]; }

    VECCORE_ATT_HOST_DEVICE
    void SetType(int type) { fType = type; }
    VECCORE_ATT_HOST_DEVICE
    void SetBr(double br) { fBr = br; }
    VECCORE_ATT_HOST_DEVICE
    void AddDaughter(int daughter) { fDaughters.push_back(daughter); }

#ifndef VECCORE_CUDA
    friend std::ostream &operator<<(std::ostream &os, const Decay &dec);
#else
    VECCORE_ATT_HOST_DEVICE
    Decay operator=(const Decay &dec) { return dec; }
#endif
  private:
    char fType;
    float fBr;
    VectorInt_t fDaughters;
  };

private:
#ifndef VECCORE_CUDA
  static void GetPart(const std::string &line, int &count, std::string &name, int &pdg, bool &matter, int &pcode,
                      std::string &pclass, int &charge, double &mass, double &width, int &isospin, int &iso3,
                      int &strange, int &flavor, int &track, int &ndecay, int &ipart, int &acode);

  static void GetDecay(const std::string &line, int &dcount, Decay &decay);
#endif
  char fName[256];          // Name
  int fPDG;                 // PDG code
  bool fMatter;             // False if antiparticle
  char fClass[256];         // Particle class
  int fPcode;               // Particle code
  float fCharge;            // Charge
  float fMass;              // Mass in GeV
  float fWidth;             // Width in GeV
  float fLife;              // Lifetime in seconds
  char fIsospin;            // Isospin
  char fIso3;               // Isospin 3
  char fStrange;            // Strangeness
  char fFlavor;             // Flavor code (?)
  char fTrack;              // Track code (?)
  unsigned char fNdecay;    // Number of decay channels
  short fCode;              // Particle code for a given MC
  VectorDecay_t fDecayList; // Decay channels
#ifndef VECCORE_CUDA
  static Map_t *fParticles; // Particle list indexed by PDG code
#endif
};
}
}
#endif
