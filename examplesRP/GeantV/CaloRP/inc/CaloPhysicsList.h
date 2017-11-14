
#ifndef CaloPhysicsList_H
#define CaloPhysicsList_H

#include "PhysicsList.h"
// for the MSCSteppingAlgorithm enums
#include "MSCModel.h"

#include <string>


namespace userapplication {

/**
 * @brief User physics list for CaloApp.
 *
 * The physics list contains the available GeantV standard EM interactions. The multiple Coulomb scattering process
 * stepping algorithm type are configurable from input arguments.
 *
 * @class   CaloPhysicsList
 * @author  M Novak
 * @date    July 2017
 */

class CaloPhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  CaloPhysicsList(const std::string &name);
  // DTR
 ~CaloPhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

  // public method to allow multiple scattering step limit configuration
  void SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping);

private:
  geantphysics::MSCSteppingAlgorithm  fMSCSteppingAlgorithm;
};

}      //  namespace userapplication


#endif // CaloPhysicsList_H
