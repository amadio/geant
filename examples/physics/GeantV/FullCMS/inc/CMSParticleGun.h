
#ifndef CMSPARTICLEGUN_H
#define CMSPARTICLEGUN_H

// from GeantV
#include "Geant/PrimaryGenerator.h"

// for GEANT_IMPL_NAMESPACE
#include "Geant/Typedefs.h"

namespace GEANT_IMPL_NAMESPACE {
  namespace geant {
    class Track;
    class TaskData;
    class EventInfo;
  }
}

// geantphysics classdef
namespace geantphysics {
  class Particle;
}


#include <string>
#include <map>

namespace cmsapp {


class CMSParticleGun : public geant::PrimaryGenerator {
public:
  // CTR DTR
  CMSParticleGun();
 ~CMSParticleGun();

  // interface methods
  virtual void InitPrimaryGenerator() {/*nothing to do*/}
  virtual geant::EventInfo NextEvent(geant::TaskData* td);
  virtual void GetTrack(int n, geant::Track& gtrack, geant::TaskData* td);

  // public setters/getters
  void SetNumPrimaryPerEvt(const int pperevt);
  int  GetNumPrimaryPerEvt() const { return fNumPrimaryPerEvt; }
  void SetPrimaryEnergy(const double ekin);
  void SetPrimaryName(const std::string& pname);
  void SetPrimaryDirection(const double pdir[3]);

  static void  Print();
  static int   GetNumberOfPrimaryTypes()  {  return gNumberCandidateParticles; }
  static int   GetPrimaryTypeIndex(const std::string& pname);
  static int   GetMaxNumberOfPrimariesPerEvent() { return gSupNumPrimaryPerEvt; }
  static const std::string&  GetPrimaryName(const int primtypeindx) { return gNameParticlesVector[primtypeindx]; }


private:
 //CMSParticleGun() = delete;
 CMSParticleGun(const CMSParticleGun&) = delete;
 CMSParticleGun &operator=(const CMSParticleGun&) = delete;

private:
  static const int          gNumberCandidateParticles;
  static const std::string  gNameParticlesVector[];
  static const std::map<std::string,int> gPrimaryNameToIndexMap;
  //
  static const int          gInfNumPrimaryPerEvt;
  static const int          gSupNumPrimaryPerEvt;
  //
  static const double       gInfBeamEnergy;
  static const double       gSupBeamEnergy;
  // these static variables stores the gun configuration just for the Print()
  static int                gNumPrimaryPerEvt;
  static double             gPrimaryEnergy;
  static std::string        gPrimaryType;
  static double             gPrimaryDir[3];
  //
  bool                      fIsUserNumPrimaryPerEvt;
  bool                      fIsUserPrimaryType;
  bool                      fIsUserPrimaryDir;
  bool                      fIsUserPrimaryEnergy;
  //
  std::string               fPrimaryParticleName;
  int                       fNumPrimaryPerEvt;
  //
  double                    fPrimaryParticleEnergy;
  //
  double                    fXPos;  // (x,y,z) position of the primary particles in internal [length] unit
  double                    fYPos;
  double                    fZPos;
  double                    fXDir;  // direction vector of the primary particles
  double                    fYDir;
  double                    fZDir;
};


}       // namespace cmsapp

#endif  // CMSPARTICLEGUN
