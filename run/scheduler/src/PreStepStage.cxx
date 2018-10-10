#include "Geant/PreStepStage.h"

#include "Geant/RunManager.h"
#include "Geant/PreStepHandler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepStage::PreStepStage(Propagator *prop) : SimulationStage(kPreStepStage, prop) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepStage::PreStepStage(const PreStepStage &other) : SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
PreStepStage &PreStepStage::operator=(const PreStepStage &other)
{
  SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage *PreStepStage::Clone() const
{
  PreStepStage *stage = new PreStepStage(*this);
  return stage;
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
Handler *PreStepStage::Select(Track *, TaskData *)
{
  // All tracks selected
  return fHandlers[0];
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
