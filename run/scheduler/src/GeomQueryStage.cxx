#include "Geant/GeomQueryStage.h"

#include "Geant/VBconnector.h"
#include "Geant/RunManager.h"
#include "Geant/GeomQueryHandler.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryStage::GeomQueryStage(Propagator *prop) : SimulationStage(kGeometryStepStage, prop) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryStage::GeomQueryStage(const GeomQueryStage &other) : SimulationStage(other) {}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryStage &GeomQueryStage::operator=(const GeomQueryStage &other)
{
  SimulationStage::operator=(other);
  return *this;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
SimulationStage *GeomQueryStage::Clone() const
{
  GeomQueryStage *stage = new GeomQueryStage(*this);
  return stage;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeomQueryStage::CreateHandlers()
{
  // Create all volume handlers.
  vector_t<Volume_t const *> &volumes = fPropagator->fRunMgr->GetVolumes();
  int nvolumes                        = fPropagator->fRunMgr->GetNvolumes();
  int threshold                       = fPropagator->fConfig->fNperBasket;
  Volume_t *vol;
  for (auto ivol = 0; ivol < nvolumes; ++ivol) {
    vol                       = (Volume_t *)volumes[ivol];
    GeomQueryHandler *handler = new GeomQueryHandler(vol, threshold, fPropagator, ivol);
    handler->SetMayBasketize(true);
    AddHandler(handler);
    assert(handler == fHandlers[ivol]);
  }
  return nvolumes;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeomQueryStage::ActivateBasketizing(bool flag)
{
  // Activate basketizers for all volumes.
  for (auto i = 0; i < GetNhandlers(); ++i)
    fHandlers[i]->ActivateBasketizing(flag);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *GeomQueryStage::Select(Track *track, TaskData *)
{
  // Retrieve the appropriate handler depending on the current volume.
  // Tracks that are killed or exit the setup should be filtered out by the relocator stage
  // and never reach this point.
  Volume_t *vol     = const_cast<Volume_t *>(track->GetVolume());
  VBconnector *link = reinterpret_cast<VBconnector *>(vol->GetBasketManagerPtr());
  return (fHandlers[link->index]);
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
