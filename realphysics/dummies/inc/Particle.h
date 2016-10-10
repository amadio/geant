#ifndef PARTICLE_H
#define PARTICLE_H

#include <string>
#include <vector>

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
  double             GetPDGCharge()    const { return fPDGCharge; }
  double             GetPDGMass()      const { return fPDGMass; }


  void ClearPhysicsProcessVector() { fPhysicsProcessVector.clear(); }
  std::vector<PhysicsProcess*>&  GetPhysicsProcessVector() {return fPhysicsProcessVector;}
  std::vector<PhysicsManagerPerParticle*>& GetPhysicsManagerPerParticleVector() {return fPMPParticle;}
  PhysicsManagerPerParticle* GetPhysicsManagerPerParticlePerRegion(int regionindx) const  {return fPMPParticle[regionindx];}


  static const Particle* GetParticleByInteralCode(unsigned int intercode) {
    if (intercode<gInternalParticleCodes.size())
      return gInternalParticleCodes[intercode];
    return nullptr;
  }

  static const std::vector<Particle*>& GetTheParticleTable() { return gTheParticleTable;}


private:
  std::string fName;
  int     fIndex;        // in the global particle table
  int     fInternalCode;
  int     fPDGCode;
  double  fPDGCharge;
  double  fPDGMass;

  std::vector<PhysicsProcess*>  fPhysicsProcessVector;  // only one and used only as temporary storage
  std::vector<PhysicsManagerPerParticle*> fPMPParticle;  // as many as regions but those having no any active processes will be nullptr

  // the particle table
  static std::vector<Particle*> gTheParticleTable;
  static std::vector<Particle*> gInternalParticleCodes;
};

} // namespace geantphysics

#endif  // PARTICLE_H
