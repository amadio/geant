
#ifndef MSCPROCESS_H
#define MSCPROCESS_H

#include "EMPhysicsProcess.h"

// from geantV
#include "Geant/Config.h"
namespace Geant {
  inline namespace GEANT_IMPL_NAMESPACE {
  class GeantTrack;
  class GeantTaskData;
}
}

namespace geantphysics {

class MSCProcess : public EMPhysicsProcess {
public:
  MSCProcess(const std::string &name);
  virtual ~MSCProcess();

  virtual void  AlongStepLimitationLength(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const ;
  virtual void  AlongStepDoIt(Geant::GeantTrack *gtrack, Geant::GeantTaskData *td) const ;

  double GetGeomMinLimit() const     { return fGeomMinLimit; }
  void   GetGeomMinLimit(double val) { fGeomMinLimit = val;  }

private:
  double fGeomMinLimit; // if the true step length is below this => no msc
};

}       // namespace geantphysics

#endif  // MSCPROCESS_H
