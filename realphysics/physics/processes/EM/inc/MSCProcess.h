
#ifndef MSCPROCESS_H
#define MSCPROCESS_H

#include "EMPhysicsProcess.h"

// from geantV
#include "Geant/Config.h"
#include "GeantTrack.h"

namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantTaskData;
}
}

namespace geantphysics {

/**
 * @brief   Base class for multiple Coulomb scattering processes.
 * @class   MSCProcess
 * @author  M Novak
 * @date    June 2017
 */


class MSCProcess : public EMPhysicsProcess {
public:
  MSCProcess(const std::string &name);
  virtual ~MSCProcess();

  virtual void  Initialize();

  virtual double  AlongStepLimitationLength(Geant::GeantTrack* /*gtrack*/, Geant::GeantTaskData* /*td*/) const ;
  // just to suppress warning since msc will implement only the second AlongStepDoIt
  using           EMPhysicsProcess::AlongStepDoIt;
  virtual void    AlongStepDoIt(Geant::GeantTrack* /*gtrack*/, Geant::GeantTaskData* /*td*/) const ;

  double GetGeomMinLimit() const     { return fGeomMinLimit; }
  void   SetGeomMinLimit(double val) { fGeomMinLimit = val;  fGeomMinLimit2 = val*val;}

private:
  double fGeomMinLimit;  // if the true step length is below this => no msc
  double fGeomMinLimit2; // square of the above
  Geant::TrackToken fMSCdata;   // Token for MSCData

};

}       // namespace geantphysics

#endif  // MSCPROCESS_H
