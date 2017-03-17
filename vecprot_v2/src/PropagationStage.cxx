#include "PropagationStage.h"

#include "GeantRunManager.h"
#include "LinearPropagationHandler.h"
#include "FieldPropagationHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PropagationStage::PropagationStage(GeantPropagator *prop)
  : SimulationStage(kPropagationStage, prop)
{
  const double bmag = prop->fConfig->fBmag;
  fHasField = bmag > 1.e-10;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int PropagationStage::CreateHandlers()
{
// Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new LinearPropagationHandler(threshold, fPropagator));
  AddHandler(new FieldPropagationHandler(threshold, fPropagator));
  
  return 2;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *PropagationStage::Select(GeantTrack *track, GeantTaskData *)
{
// Retrieve the appropriate handler depending on the track charge
  if (!fHasField || track->fCharge == 0)
    return fHandlers[0];
  return fHandlers[1];
}

} // GEANT_IMPL_NAMESPACE
} // Geant
