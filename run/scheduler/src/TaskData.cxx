#include "Geant/TaskData.h"
#include "Geant/Basket.h"
#include "Geant/BasketCounters.h"
#include "StackLikeBuffer.h"
#include "Geant/Propagator.h"
#include "Geant/TrackManager.h"
#include "TrackGeo.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "Geant/SimulationStage.h"
#include "TrackStat.h"

#include "Geant/GUFieldPropagator.h"
#include "Geant/VVectorField.h"

#ifdef USE_ROOT
#include "TRandom.h"
#endif

using std::min;

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
TaskData::TaskData(size_t nthreads, int maxPerBasket)
{
  // Constructor
  fNthreads = nthreads;
  fSizeBool = fSizeDbl = fSizeInt = 5 * maxPerBasket;
  fBoolArray                      = new bool[fSizeBool];
  fDblArray                       = new double[fSizeDbl];
  fIntArray                       = new int[fSizeInt];
  fPath                           = VolumePath_t::MakeInstance(TrackDataMgr::GetInstance()->GetMaxDepth());
  fPathV                          = new VolumePath_t *[4 * maxPerBasket];
  fNextpathV                      = new VolumePath_t *[4 * maxPerBasket];
  fTrack                          = Track::MakeInstance();
  fGeoTrack                       = new TrackGeo_v(4 * maxPerBasket);
  fRndm                           = new RngWrapper;
  fQshare                         = new queue_t(1 << 13);
}

//______________________________________________________________________________
TaskData::~TaskData()
{
  // Destructor
  Track::ReleaseInstance(fTrack);
  delete[] fBoolArray;
  delete[] fDblArray;
  delete[] fIntArray;
  delete[] fPathV;
  delete[] fNextpathV;
  delete fRndm;
  VolumePath_t::ReleaseInstance(fPath);
  delete fStackBuffer;
  delete fShuttleBasket;
  delete fBvector;
  delete fStat;
  for (auto basket : fStageBuffers)
    delete basket;
  fStageBuffers.clear();
  for (auto stage : fStages) {
    if (stage->HasLocalHandlers()) {
      stage->DeleteLocalHandlers();
      delete stage;
    }
  }
  fStages.clear();
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TaskData::AttachPropagator(Propagator *prop, int node)
{
  // Attach to a given propagator and a given NUMA node.
  if (fPropagator) {
    assert(fPropagator == prop);
    fNode = node;
    return;
  }
  fPropagator    = prop;
  fNode          = node;
  bool usenuma   = prop->fConfig->fUseNuma;
  fShuttleBasket = usenuma ? new Basket(1000, 0, node) : new Basket(1000, 0);
  fBvector       = usenuma ? new Basket(256, 0, node) : new Basket(256, 0);
  for (int i = 0; i <= int(kSteppingActionsStage); ++i)
    fStageBuffers.push_back(usenuma ? new Basket(1000, 0, node) : new Basket(1000, 0));
  fStackBuffer = new StackLikeBuffer(prop->fConfig->fNstackLanes, this);
  fStackBuffer->SetStageBuffer(fStageBuffers[0]);
  fBlock = fPropagator->fTrackMgr->GetNewBlock();
  for (size_t stage = 0; stage < kNstages; ++stage) {
    // Check if the stage has thread local handlers
    if (!prop->fStages[stage]->IsBasketized() || !prop->fStages[stage]->HasLocalHandlers()) {
      fStages.push_back(prop->fStages[stage]);
    } else {
      // Replace the stage with a thread local one
      // Printf("...cloning stage %s into task data %d", prop->fStages[stage]->GetName(), fTid);
      auto clone_stage = prop->fStages[stage]->Clone();
      clone_stage->ReplaceLocalHandlers();
      fStages.push_back(clone_stage);
    }
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TaskData::CreateStageCounters(Propagator *prop)
{
  for (size_t stage = 0; stage < kNstages; ++stage) {
    fCounters[stage] = new BasketCounters(prop->fStages[stage]->GetNhandlers());
  }
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track &TaskData::GetNewTrack()
{
  size_t index;
  if (fBlock->IsDistributed()) {
    fBlock = fPropagator->fTrackMgr->GetNewBlock();
    // printf("== New block: %d (%d) current=%d used=%d\n",
    //       fBlock->GetId(), fBlock->GetNode(), fBlock->GetCurrent(), fBlock->GetUsed());
    assert(fBlock->GetCurrent() == 0 && fBlock->GetUsed() == 0);
  }
  Track *track = fBlock->GetObject(index);
  track->Reset(*fTrack);
  track->SetBindex(index);
  return *track;
  //  Track &track = fPropagator->fTrackMgr->GetTrack();
  //  index = track.BIndex();
  //  track.Reset(*fTrack);
  //  track.SetBindex(index);
  //  return track;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void TaskData::ReleaseTrack(Track &track)
{
  fPropagator->fTrackMgr->ReleaseTrack(track);
}

//______________________________________________________________________________
void TaskData::InspectStages(int istage) const
{
  Printf("** Thread %d: **", fTid);
  for (auto stage : fPropagator->fStages) {
    if (stage->GetId() == istage)
      Printf("*** -> %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());
    else
      Printf("***    %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());
  }
}

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
