#include "GeomQueryStage.h"

#include "GeantRunManager.h"
#include "GeomLengthQuery.h"

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
int GeomQueryStage::CreateSelectors()
{
// Create all volume selectors.
  vector_t<Volume_t const *> volumes = fPropagator->fRunMgr->GetVolumes();
  int nvolumes = fPropagator->fRunMgr->GetNvolumes();
  int threshold = fPropagator->fConfig->fNperBasket;
  Volume_t *vol;
  for (auto ivol = 0; ivol < nvolumes; ++ivol) {
    vol = (Volume_t *)volumes[ivol];
    GeomLengthQuery *selector = new GeomLengthQuery(vol, threshold, fPropagator, ivol);
    AddSelector(selector);
    assert(selector == fSelectors[ivol]);
  }
  return nvolumes;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Selector *GeomQueryStage::Select(GeantTrack *track)
{
// Retrieve the appropriate selector depending on the current volume.
// Tracks that are killed or exit the setup should be filtered out by the relocator stage
// and never reach this point.
#ifdef USE_VECGEOM_NAVIGATOR
  Volume_t *vol = const_cast<Volume_t *>(track->fPath->Top()->GetLogicalVolume());
  VBconnector *link = reinterpret_cast<VBconnector *>(vol->GetBasketManagerPtr());
#else
  Volume_t *vol = track->fPath->GetCurrentNode()->GetVolume();
  VBconnector *link = reinterpret_cast<VBconnector *>(vol->GetFWExtension());
#endif
  return ( fSelectors[link->index] );
}

} // GEANT_IMPL_NAMESPACE
} // Geant
