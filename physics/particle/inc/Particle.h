#ifndef Particle_H
#define Particle_H

#include "base/Global.h"
#include "Geant/Config.h"

#ifdef GEANT_NVCC
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

#ifdef GEANT_NVCC
class Particle;
extern GEANT_CUDA_DEVICE_CODE vecgeom::map<int, Particle> *fParticlesDev; // Particle list indexed by PDG code
extern vecgeom::map<int, Particle> *fParticlesHost;                           // Particle list indexed by PDG code
#endif

class Particle {
public:
  class Decay;
#ifdef GEANT_NVCC
  using Map_t         = vecgeom::map<int, Particle>;
  using VectorDecay_t = vecgeom::Vector<Decay>;
  using VectorInt_t   = vecgeom::Vector<int>;
#else
  using Map_t         = std::map<int, Particle>;
  using VectorDecay_t = std::vector<Decay>;
  using VectorInt_t   = std::vector<int>;
#endif

  GEANT_CUDA_BOTH_CODE
  Particle();
  GEANT_CUDA_BOTH_CODE
  Particle(const char *name, int pdg, bool matter, const char *pclass, int pcode, double charge, double mass,
           double width, int isospin, int iso3, int strange, int flavor, int track, int code = -1);

  GEANT_CUDA_BOTH_CODE
  Particle(const Particle &other);

  GEANT_CUDA_BOTH_CODE
  Particle &operator=(const Particle &part);

  GEANT_CUDA_BOTH_CODE
  static void CreateParticles();
  GEANT_CUDA_BOTH_CODE
  const char *Name() const { return fName; }
  GEANT_CUDA_BOTH_CODE
  int PDG() const { return fPDG; }
  GEANT_CUDA_BOTH_CODE
  bool Matter() const { return fMatter; }
  GEANT_CUDA_BOTH_CODE
  double Mass() const { return fMass; }
  const char *Class() const { return fClass; }
  GEANT_CUDA_BOTH_CODE
  int Pcode() const { return fPcode; }
  GEANT_CUDA_BOTH_CODE
  double Charge() const { return fCharge; }
  GEANT_CUDA_BOTH_CODE
  double Width() const { return fWidth; }
  GEANT_CUDA_BOTH_CODE
  int Isospin() const { return fIsospin; }
  GEANT_CUDA_BOTH_CODE
  int Iso3() const { return fIso3; }
  GEANT_CUDA_BOTH_CODE
  int Strange() const { return fStrange; }
  GEANT_CUDA_BOTH_CODE
  int Flavor() const { return fFlavor; }
  GEANT_CUDA_BOTH_CODE
  int Track() const { return fTrack; }
  GEANT_CUDA_BOTH_CODE
  int Ndecay() const { return fNdecay; }
  GEANT_CUDA_BOTH_CODE
  int Code() const { return fCode; }

  GEANT_CUDA_BOTH_CODE
  void SetCode(int code) { fCode = code; }

  const VectorDecay_t &DecayList() const { return fDecayList; }
#ifndef GEANT_NVCC
  static void ReadFile(std::string infilename, bool outfile=false);
#endif

#ifndef GEANT_NVCC
  static const Particle &GetParticle(int pdg)
  {
    if (fParticles->find(pdg) != fParticles->end()) return (*fParticles)[pdg];
    static Particle p;
    std::cout << __func__ << "::pdg:" << pdg << " does not exist" << std::endl;
    return p;
  }
#else
  static const Particle &GetParticle(int pdg)
  {
    if (fParticlesHost->find(pdg) != fParticlesHost->end()) return (*fParticlesHost)[pdg];
    static Particle p;
    printf(" pdg %d does not exist\n", pdg);
    return p;
  }
  GEANT_CUDA_DEVICE_CODE
  static const Particle &GetParticleDev(int pdg)
  {
    if (fParticlesDev->find(pdg) != fParticlesDev->end()) return (*fParticlesDev)[pdg];
    // Particle p;
    printf(" pdg %d does not exist\n", pdg);
    return (*fParticlesDev)[1];
  }
#endif

#ifndef GEANT_NVCC
  void NormDecay();

  friend std::ostream &operator<<(std::ostream &os, const Particle &part);

#endif
  GEANT_CUDA_BOTH_CODE
  void AddDecay(const Decay &decay)
  {
    fDecayList.push_back(decay);
    fNdecay = fDecayList.size();
  }
  GEANT_CUDA_BOTH_CODE
  static const Map_t &Particles()
  {
#ifndef GEANT_NVCC
    return *fParticles;
#else
#ifndef GEANT_CUDA_DEVICE_CODE
    return *fParticlesHost;
#else
    return *fParticlesDev;
#endif
#endif
  }

  class Decay {
  public:
    GEANT_CUDA_BOTH_CODE
    Decay() : fType(0), fBr(0) {}
    GEANT_CUDA_BOTH_CODE
    Decay(const Decay &other) : fType(other.fType), fBr(other.fBr), fDaughters(other.fDaughters) {}
    GEANT_CUDA_BOTH_CODE
    Decay(int type, double br, const VectorInt_t &daughters) : fType(type), fBr(br), fDaughters(daughters) {}
    GEANT_CUDA_BOTH_CODE
    void Clear()
    {
      fType = 0;
      fBr   = 0;
      fDaughters.clear();
    }

    GEANT_CUDA_BOTH_CODE
    int Type() const { return fType; }
    GEANT_CUDA_BOTH_CODE
    double Br() const { return fBr; }
    GEANT_CUDA_BOTH_CODE
    const VectorInt_t &Daughters() const { return fDaughters; }
    GEANT_CUDA_BOTH_CODE
    int NDaughters() const { return fDaughters.size(); }
    GEANT_CUDA_BOTH_CODE
    int Daughter(int i) const { return fDaughters[i]; }

    GEANT_CUDA_BOTH_CODE
    void SetType(int type) { fType = type; }
    GEANT_CUDA_BOTH_CODE
    void SetBr(double br) { fBr = br; }
    GEANT_CUDA_BOTH_CODE
    void AddDaughter(int daughter) { fDaughters.push_back(daughter); }

#ifndef GEANT_NVCC
    friend std::ostream &operator<<(std::ostream &os, const Decay &dec);
#else
    GEANT_CUDA_BOTH_CODE
    Decay operator=(const Decay &dec) { return dec; }
#endif
  private:
    char fType;
    float fBr;
    VectorInt_t fDaughters;
  };

private:
#ifndef GEANT_NVCC
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
#ifndef GEANT_NVCC
  static Map_t *fParticles; // Particle list indexed by PDG code
#endif
};
}
}
#endif
