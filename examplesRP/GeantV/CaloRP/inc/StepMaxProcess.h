
#ifndef STEPMAXPROCESS_H
#define STEPMAXPROCESS_H

#include "PhysicsProcess.h"

namespace userapplication {

/**
 * @brief Simple user defined physics-process to limit the (true) step length of e-/e+.
 *
 * The process implements only the PostStepLimitationLength (since it will act as discrete process) physics process
 * interface method. It will make sure that a particle, to which this process is assigned to in the physics list, never
 * makes longer (true) step length than a given value. This value can be set by using the SetMaxStep() method. When this
 * step-max process limits the step, it will act as a delta process: nothing happens.
 *
 * @class   StepMaxProcess
 * @author  M Novak
 * @date    July 2017
 */

class StepMaxProcess : public geantphysics::PhysicsProcess {
public:
  StepMaxProcess(const std::string &name = "Step-max-process");
 ~StepMaxProcess();

  // interface method
 virtual double PostStepLimitationLength(Geant::GeantTrack * /*track*/, Geant::GeantTaskData * /*td*/, bool haseloss=false);


 void SetMaxStep(double val) { fMaxStep = val; }

private:
  double  fMaxStep;

};

}        // namespace userapplication

#endif   // STEPMAXPROCESS_H
