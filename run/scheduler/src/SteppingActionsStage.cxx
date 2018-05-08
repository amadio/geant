#include "Geant/SteppingActionsStage.h"

#include "Geant/RunManager.h"
#include "Geant/SteppingActionsHandler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsStage::SteppingActionsStage(Propagator *prop) : SimulationStage(kSteppingActionsStage, prop) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsStage::SteppingActionsStage(const SteppingActionsStage &other) : SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SteppingActionsStage &SteppingActionsStage::operator=(const SteppingActionsStage &other)
{
  SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage *SteppingActionsStage::Clone() const
{
  SteppingActionsStage *stage = new SteppingActionsStage(*this);
  return stage;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int SteppingActionsStage::CreateHandlers()
{
  // Create all volume handlers.
  int threshold = fPropagator->fConfig->fNperBasket;
  bool verbose  = fPropagator->fConfig->fSteppingVerbose;
  auto stepHndl = new SteppingActionsHandler(threshold, fPropagator);
  stepHndl->SetVerbosity(verbose);
  AddHandler(stepHndl);
  if( verbose ) 
    std::cout << "SteppingACtionsStage: set verbose = " << ( verbose ? " ON" : " OFF" ) << std::endl;

  return 1;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *SteppingActionsStage::Select(Track *, TaskData *)
{
  // Retrieve the appropriate handler depending on the track charge
  return fHandlers[0];
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
