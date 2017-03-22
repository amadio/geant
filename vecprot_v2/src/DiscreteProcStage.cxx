#include "DiscreteProcStage.h"

#include "GeantTaskData.h"
#include "PhysicsProcessOld.h"
#include "DiscreteProcHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
DiscreteProcStage::DiscreteProcStage(GeantPropagator *prop)
  : SimulationStage(kDiscreteProcStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int DiscreteProcStage::CreateHandlers()
{
// Create the discrete process handlers
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new DiscreteProcHandler(threshold, fPropagator));
  
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *DiscreteProcStage::Select(GeantTrack *track, GeantTaskData *td)
{
// Select tracks with limiting discrete process
  if (track->fStatus == kPhysics && track->fEindex == 1000) {
    // reset number of interaction length left
    track->fNintLen = -1;
    // Invoke PostStepTypeOfIntrActSampling
    fPropagator->Process()->PostStepTypeOfIntrActSampling(track, td);

    return ( GetHandler(track->fProcess) );
  }
  return nullptr;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
