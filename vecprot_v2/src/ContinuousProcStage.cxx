#include "ContinuousProcStage.h"

#include "ContinuousProcHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
ContinuousProcStage::ContinuousProcStage(GeantPropagator *prop)
  : SimulationStage(kContinuousProcStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int ContinuousProcStage::CreateHandlers()
{
// Create the continuous process handlers
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new ContinuousProcHandler(threshold, fPropagator));
  
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *ContinuousProcStage::Select(GeantTrack *track, GeantTaskData *)
{
// Retrieve the appropriate handler depending on the track charge
  if (track->Charge() == 0)
    return nullptr;
  return fHandlers[0];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
