#include "PreStepStage.h"

#include "GeantRunManager.h"
#include "PreStepHandler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepStage::PreStepStage(GeantPropagator *prop)
  : SimulationStage(kPreStepStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int PreStepStage::CreateHandlers()
{
// Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new PreStepHandler(threshold, fPropagator));
  
  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *PreStepStage::Select(GeantTrack *, GeantTaskData *)
{
// All tracks selected
  return fHandlers[0];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
