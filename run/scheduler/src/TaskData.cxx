#include "TaskData.h"
#include "Basket.h"
#include "BasketCounters.h"
#include "StackLikeBuffer.h"
#include "Propagator.h"
#include "TrackManager.h"
#include "TrackGeo.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "SimulationStage.h"
#include "TrackStat.h"

#include "GUFieldPropagator.h"
#include "VVectorField.h"

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
  fBoolArray = new bool[fSizeBool];
  fDblArray = new double[fSizeDbl];
  fIntArray = new int[fSizeInt];
  fPath = VolumePath_t::MakeInstance(TrackDataMgr::GetInstance()->GetMaxDepth());
  fPathV = new VolumePath_t*[4*maxPerBasket];
  fNextpathV = new VolumePath_t*[4*maxPerBasket];
  fTrack = Track::MakeInstance();
  fGeoTrack = new TrackGeo_v(4*maxPerBasket);
  fRndm = new vecgeom::RNG; // what about the seed?
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
TaskData::TaskData(void *addr, size_t nthreads, int maxPerBasket, Propagator *prop)
    : fPropagator(prop), fNthreads(nthreads)
{
  // Constructor
  char *buffer = (char*)addr;
  buffer += Track::round_up_align(sizeof(TaskData));
  buffer = Track::round_up_align(buffer);

  fPath = VolumePath_t::MakeInstanceAt(TrackDataMgr::GetInstance()->GetMaxDepth(), (void*)buffer);
  fPathV = new VolumePath_t*[4*maxPerBasket];
  fNextpathV = new VolumePath_t*[4*maxPerBasket];
  fGeoTrack = TrackGeo_v::MakeInstanceAt(buffer, 4*maxPerBasket);
  buffer += TrackGeo_v::SizeOfInstance(4*maxPerBasket);
  buffer += VolumePath_t::SizeOfInstance(TrackDataMgr::GetInstance()->GetMaxDepth());
  buffer = Track::round_up_align(buffer);

  // Previous, the size was hard coded to 1024, '4' is a guess on the max number
  // of produced particles ...
  fSizeInt = fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new (buffer) bool[fSizeBool];
  buffer += fSizeBool*sizeof(bool);
  fDblArray = new (buffer) double[fSizeDbl];
  buffer += fSizeDbl*sizeof(double);
  fIntArray = new (buffer) int[fSizeInt];
  buffer += fSizeInt*sizeof(int);

  fTrack = Track::MakeInstance();

  fRndm = &vecgeom::RNG::Instance();
}

//______________________________________________________________________________
TaskData::~TaskData()
{
// Destructor
  Track::ReleaseInstance(fTrack);
  delete[] fBoolArray;
  delete[] fDblArray;
  delete[] fIntArray;
  delete [] fPathV;
  delete [] fNextpathV;
  delete fRndm;
  VolumePath_t::ReleaseInstance(fPath);
  delete fStackBuffer;
  delete fShuttleBasket;
  delete fBvector;
  delete fStat;
  for (auto basket : fStageBuffers)
    delete basket;
  fStageBuffers.clear();
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
  fPropagator = prop;
  fNode = node;
  bool usenuma = prop->fConfig->fUseNuma;
  fShuttleBasket = usenuma ? new Basket(1000, 0, node) : new Basket(1000, 0);
  fBvector = usenuma ? new Basket(256, 0, node) : new Basket(256, 0);
  for (int i=0; i<=int(kSteppingActionsStage); ++i)
    fStageBuffers.push_back(usenuma ? new Basket(1000, 0, node) : new Basket(1000, 0));
  fStackBuffer = new StackLikeBuffer(prop->fConfig->fNstackLanes, this);
  fStackBuffer->SetStageBuffer(fStageBuffers[0]);
  fBlock = fPropagator->fTrackMgr->GetNewBlock();
  for (size_t stage = 0; stage < kNstages; ++stage)
    fCounters[stage] = new BasketCounters(prop->fStages[stage]->GetNhandlers());
  // std::cerr<<"AttachProp(): prop="<< prop <<", node="<< node <<", fBlock="<< fBlock <<"\n";
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
TaskData *TaskData::MakeInstanceAt(void *addr, size_t nTracks, int maxPerBasket, Propagator *prop)
{
   // Track MakeInstance based on a provided single buffer.
   return new (addr) TaskData(addr, nTracks, maxPerBasket, prop);
}


//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t TaskData::SizeOfInstance(size_t /*nthreads*/, int maxPerBasket)
{
   // @brief return the contiguous memory size needed to hold a TaskData

   const size_t bufSize = 5; // See constructor!

   size_t need = sizeof(TaskData) // vecgeom::DevicePtr<geant::cuda::TaskData>::SizeOf()
      + Track::round_up_align(bufSize*maxPerBasket*(sizeof(bool)+sizeof(double)+sizeof(int)))
      + Track::round_up_align(VolumePath_t::SizeOfInstance(TrackDataMgr::GetInstance()->GetMaxDepth()));
   return Track::round_up_align(need);
}


#ifndef VECCORE_CUDA

//______________________________________________________________________________
void TaskData::RecycleBasket(Basket *b)
{
  // Recycle a basket.
  fBPool.push_back(b);
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
Track &TaskData::GetNewTrack()
{
  size_t index;
  if (fBlock->IsDistributed()) {
    fBlock = fPropagator->fTrackMgr->GetNewBlock();
    //printf("== New block: %d (%d) current=%d used=%d\n",
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
void TaskData::ReleaseTrack(Track &track) {
  fPropagator->fTrackMgr->ReleaseTrack(track);
}

//______________________________________________________________________________
void TaskData::InspectStages(int istage)
{
  Printf("** Thread %d: **", fTid);
  for (auto stage : fPropagator->fStages) {
    if (stage->GetId() == istage)
      Printf("*** -> %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());
    else
      Printf("***    %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());    
  }
}


#endif

} // GEANT_IMPL_NAMESPACE
} // geant
