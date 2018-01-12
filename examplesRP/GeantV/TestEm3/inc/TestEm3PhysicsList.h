
#ifndef TESTEM3PHYSICSLIST_H
#define TESTEM3PHYSICSLIST_H

#include "PhysicsList.h"
// for the MSCSteppingAlgorithm enums
#include "MSCModel.h"

#include <string>


namespace userapplication {

/**
 * @brief User physics list for TestEm3App.
 *
 * The physics list contains the available GeantV standard EM interactions. The multiple Coulomb scattering process
 * stepping algorithm type are configurable from input arguments.
 *
 * @class   TestEm3PhysicsList
 * @author  M Novak
 * @date    July 2017
 */

class TestEm3PhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  TestEm3PhysicsList(const std::string &name);
  // DTR
 ~TestEm3PhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

  // public method to allow multiple scattering step limit configuration
  void SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping);

private:
  geantphysics::MSCSteppingAlgorithm  fMSCSteppingAlgorithm;
};

}      //  namespace userapplication


#endif // TESTEM3PHYSICSLIST_H