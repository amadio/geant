#include "GeomQueryStage.h"

#include "GeantRunManager.h"
#include "GeomQueryHandler.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeomQueryStage::GeomQueryStage(GeantPropagator *prop)
  : SimulationStage(kGeometryStepStage, prop)
{
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
int GeomQueryStage::CreateHandlers()
{
// Create all volume handlers.
  vector_t<Volume_t const *> &volumes = fPropagator->fRunMgr->GetVolumes();
  int nvolumes = fPropagator->fRunMgr->GetNvolumes();
  int threshold = fPropagator->fConfig->fNperBasket;
  Volume_t *vol;
  for (auto ivol = 0; ivol < nvolumes; ++ivol) {
    vol = (Volume_t *)volumes[ivol];
    GeomQueryHandler *handler = new GeomQueryHandler(vol, threshold, fPropagator, ivol);
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
  for ( auto i = 0; i < GetNhandlers(); ++i)
    fHandlers[i]->ActivateBasketizing(true);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Handler *GeomQueryStage::Select(GeantTrack *track, GeantTaskData *)
{
// Retrieve the appropriate handler depending on the current volume.
// Tracks that are killed or exit the setup should be filtered out by the relocator stage
// and never reach this point.
#ifdef USE_VECGEOM_NAVIGATOR
  Volume_t *vol = const_cast<Volume_t *>(track->fPath->Top()->GetLogicalVolume());
  VBconnector *link = reinterpret_cast<VBconnector *>(vol->GetBasketManagerPtr());
#else
  Volume_t *vol = track->fPath->GetCurrentNode()->GetVolume();
  VBconnector *link = reinterpret_cast<VBconnector *>(vol->GetFWExtension());
#endif
  return ( fHandlers[link->index] );
}

} // GEANT_IMPL_NAMESPACE
} // Geant
