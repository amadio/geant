
#ifndef USERPHYSICSLIST_H
#define USERPHYSICSLIST_H

#include "PhysicsList.h"
// for the MSCSteppingAlgorithm enums
#include "MSCModel.h"

#include <string>


namespace userapplication {

/**
 * @brief User physics list for TeatEm5.
 *
 * The physics list contains the available GeantV standard EM interactions and a user defined StepMaxProcess for e-/e+.
 * The step limit value of the this StepMaxProcess and the multiple Coulomb scattering process stepping algorithm type
 * are configurable from input arguments.
 *
 * @class   UserPhysicsList
 * @author  M Novak
 * @date    July 2017
 */

class UserPhysicsList : public geantphysics::PhysicsList {
public:
  // CTR
  UserPhysicsList(const std::string &name);
  // DTR
 ~UserPhysicsList();
  // interface method to assigne physics-process to particles
  virtual void Initialize();

  // public method to allow multiple scattering step limit configuration
  void SetMSCStepLimit(geantphysics::MSCSteppingAlgorithm stepping);
  // public method to allow configuration of the user-defined step-max process
  void SetStepMaxValue(double val);

private:
  geantphysics::MSCSteppingAlgorithm  fMSCSteppingAlgorithm;
  double                              fStepMaxValue;
};

}      //  namespace userapplication


#endif // USERPHYSICSLIST_H
