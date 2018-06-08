//====================================================================================================================
// TestEm3FastSimProcess.h - Geant-V Prototype
/**
 * @file TestEm3FastSimProcess.h
 * @brief Example Fast Sim process.
 *
 * The class TestEm3FastSimProcess is an example fast sim processes.
 *
 * @author W. Pokorski (May 2018)
 */
//====================================================================================================================

#ifndef TESTEM3FASTSIM_PROCESS
#define TESTEM3FASTSIM_PROCESS

#include <string>
#include <vector>
#include "Geant/FastSimProcess.h"

namespace geantphysics {
inline namespace GEANT_IMPL_NAMESPACE {
class Isotope;
class Material;
class Element;
}
}

/**
 * @brief Class TestEm3FastSimProcess
 */

namespace userapplication {
  
  class TestEm3FastSimProcess : public geantphysics::FastSimProcess {
public:
  /** @brief TestEm3FastSimProcess default constructor */
  TestEm3FastSimProcess();

  /** @brief TestEm3FastSimProcess complete constructor */
  TestEm3FastSimProcess(const std::vector<int> &particlecodevec);

  /** @brief TestEm3FastSimProcess destructor */
  virtual ~TestEm3FastSimProcess();

  // The methods below are those inherited from PhysicsProcess

  /** Method that returns "true" ("false") if the specified GV particle code is (not) accepted by this process */
  virtual bool IsApplicable(geant::Track *track) const;

  /** Main method that calls fast sim process */
  virtual int FastSimDoIt(geantphysics::LightTrack &track, geant::TaskData *td);



};
}

#endif
