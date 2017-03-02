#include "XSecSamplingStage.h"

#include "GeantRunManager.h"
#include "XSecSamplingHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
XSecSamplingStage::XSecSamplingStage(GeantPropagator *prop)
  : SimulationStage(kXSecSamplingStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int XSecSamplingStage::CreateHandlers()
{
// Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new XSecSamplingHandler(threshold, fPropagator));
  
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *XSecSamplingStage::Select(GeantTrack *, GeantTaskData *)
{
// Retrieve the appropriate handler depending on the track charge
  return fHandlers[0];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
