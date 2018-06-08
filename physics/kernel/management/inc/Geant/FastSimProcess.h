//====================================================================================================================
// FastSimProcess.h - Geant-V Prototype
/**
 * @file FastSimProcess.h
 * @brief Base class for all fast sim processes in the Geant-V prototype.
 *
 * The class FastSimProcess is the base class for all fast sim processes.
 *
 * @author W. Pokorski (May 2018)
 */
//====================================================================================================================

#ifndef FASTSIM_PROCESS
#define FASTSIM_PROCESS

#include <string>
#include <vector>
#include "Geant/PhysicsProcess.h"

namespace geantphysics {

inline namespace GEANT_IMPL_NAMESPACE {
class Isotope;
class Material;
class Element;
}

/**
 * @brief Class FastSimProcess
 */
class FastSimProcess : public PhysicsProcess {
public:
  /** @brief FastSimProcess default constructor */
  FastSimProcess();

  FastSimProcess(const std::string &name);

  /** @brief FastSimProcess complete constructor */
  FastSimProcess(const std::string &name, const std::vector<int> &particlecodevec);

  /** @brief FastSimProcess destructor */
  virtual ~FastSimProcess();

  // The methods below are those inherited from PhysicsProcess

  /** Method that returns "true" ("false") if the specified GV particle code is (not) accepted by this process */
  virtual bool IsApplicable(geant::Track *track) const;

  /** Main method that calls fast sim process */
  virtual int FastSimDoIt(LightTrack &track, geant::TaskData *td);

  //--- Getters ---

  /** Method that returns the vector of GV particle codes for the allowed particles */
  const std::vector<int> &GetParticleCodeVec() const { return fParticleCodeVec; }

  //--- Setters ---

  /** Method that sets the GV particle codes of the allowed particles */
  void SetParticleCodeVec(const std::vector<int> &particlecodevec)
  {
    fParticleCodeVec.clear();
    for (size_t i = 0; i < particlecodevec.size(); i++) {
      fParticleCodeVec.push_back(particlecodevec[i]);
    }
  }

private:
  std::vector<int> fParticleCodeVec;         /** Vector of GV particle codes for the allowed particles */
};

} // end of namespace geant

#endif
