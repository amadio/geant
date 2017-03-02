#include "SteppingActionsStage.h"

#include "GeantRunManager.h"
#include "SteppingActionsHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsStage::SteppingActionsStage(GeantPropagator *prop)
  : SimulationStage(kSteppingActionsStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SteppingActionsStage::CreateHandlers()
{
// Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new SteppingActionsHandler(threshold, fPropagator));
  
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *SteppingActionsStage::Select(GeantTrack *, GeantTaskData *)
{
// Retrieve the appropriate handler depending on the track charge
  return fHandlers[0];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
