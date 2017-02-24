#include "DiscreteProcStage.h"

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
Handler *DiscreteProcStage::Select(GeantTrack *track)
{
// Retrieve the appropriate handler depending on the track charge
  return fHandlers[0];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
