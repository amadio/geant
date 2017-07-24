#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <vector>
#include <map>

namespace geantphysics {
// forward declarations
class PhysicsProcess;
class PhysicsManagerPerParticle;
/**
 * @brief   Base class to describe particles static properties.
 * @class   Particle
 * @author  M Novak, A Ribon
 * @date    april 2016
 */
class Particle {
public:
  Particle(const std::string &name, int pdgcode, int intcode, double mass, double charge);
 ~Particle(){}

  const std::string& GetName()         const { return fName;  }
  int                GetIndex()        const { return fIndex; }
  int                GetInternalCode() const { return fInternalCode; }
  int                GetPDGCode()      const { return fPDGCode; }
  double             GetPDGMass()      const { return fPDGMass; }
  double             GetPDGCharge()    const { return fPDGCharge; }


  void ClearPhysicsProcessVector() { fPhysicsProcessVector.clear(); }
  std::vector<PhysicsProcess*>&  GetPhysicsProcessVector() {return fPhysicsProcessVector;}
  std::vector<PhysicsManagerPerParticle*>& GetPhysicsManagerPerParticleVector() {return fPMPParticle;}
  PhysicsManagerPerParticle* GetPhysicsManagerPerParticlePerRegion(int regionindx) const  {return fPMPParticle[regionindx];}


  static const Particle* GetParticleByInternalCode(unsigned int intercode) {
    if (intercode<gInternalParticleCodes.size())
      return gInternalParticleCodes[intercode];
    return nullptr;
  }

  static const Particle* GetParticleByPDGCode(unsigned int pdgcode) {
    auto search = gPDGtoInternalCode.find(pdgcode);
    if (search != gPDGtoInternalCode.end())
      return gInternalParticleCodes[search->second];
    return nullptr;
  }

  static const Particle* GetParticleByName(const std::string pname) {
    auto search = gNametoParticle.find(pname);
    if (search != gNametoParticle.end())
      return search->second;
    return nullptr;
  }

  static const std::vector<Particle*>& GetTheParticleTable() { return gTheParticleTable;}


private:
  std::string fName;
  int     fIndex;        // in the global particle table
  int     fInternalCode;
  int     fPDGCode;
  double  fPDGMass;
  double  fPDGCharge;

  std::vector<PhysicsProcess*>  fPhysicsProcessVector;  // only one and used only as temporary storage
  std::vector<PhysicsManagerPerParticle*> fPMPParticle;  // as many as regions but those having no any active processes will be nullptr

  // the particle table
  static std::vector<Particle*> gTheParticleTable;
  static std::vector<Particle*> gInternalParticleCodes;

  // map of PDG codes to internal codes
  static std::map<unsigned int, unsigned int> gPDGtoInternalCode;
  // map of name to particle ptr
  static std::map<const std::string, Particle*> gNametoParticle;

};

} // namespace geantphysics

#endif  // PARTICLE_H
