
#ifndef POSITRONANNIHILATIONPROCESS_H
#define POSITRONANNIHILATIONPROCESS_H

#include "EMPhysicsProcess.h"

#include <string>

namespace geantphysics {

/**
 * @brief   Pre-prepared physics process to describe positron electron annihilation into 2 gamma.
 * @class   PositronAnnihilationProcess
 * @author  M Novak
 * @date    January 2018
 */


class PositronAnnihilationProcess : public EMPhysicsProcess {
public:
  PositronAnnihilationProcess(const std::string &name = "Annihilation");

  virtual void   Initialize();

  virtual double AverageLifetime(const LightTrack &track) const;

  virtual int    AtRestDoIt(LightTrack &track, geant::GeantTaskData *td);
};

}        // namespace geantphysics

#endif   // POSITRONANNIHILATIONPROCESS_H
