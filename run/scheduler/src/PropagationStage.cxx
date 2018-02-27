#include "Geant/PropagationStage.h"

#include "Geant/RunManager.h"
#include "Geant/FieldLookup.h"
#include "Geant/LinearPropagationHandler.h"
#include "Geant/FieldPropagationHandler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PropagationStage::PropagationStage(Propagator *prop) : SimulationStage(kPropagationStage, prop)
{
  auto fldConfig = FieldLookup::GetFieldConfig();
  assert(fldConfig != nullptr);
  fHasField = fldConfig->FieldExists();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int PropagationStage::CreateHandlers()
{
  // Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  AddHandler(new LinearPropagationHandler(threshold, fPropagator));
  auto handler = new FieldPropagationHandler(threshold, fPropagator);
  handler->SetMayBasketize(true);
  AddHandler(handler);

  return 2;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *PropagationStage::Select(Track *track, TaskData *)
{
  // Retrieve the appropriate handler depending on the track charge
  if (!fHasField || track->Charge() == 0) return fHandlers[0];
  return fHandlers[1];
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
