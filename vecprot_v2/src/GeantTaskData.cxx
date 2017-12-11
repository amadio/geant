#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "Basket.h"
#include "BasketCounters.h"
#include "StackLikeBuffer.h"
#include "GeantPropagator.h"
#include "TrackManager.h"
#include "GeantTrackGeo.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "SimulationStage.h"
#include "TrackStat.h"

#include "GUFieldPropagator.h"
#include "VScalarField.h"

#ifdef USE_ROOT
#include "TRandom.h"
#endif

using std::min;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantTaskData::GeantTaskData(size_t nthreads, int maxPerBasket)
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
  fTrack = GeantTrack::MakeInstance();
  fGeoTrack = new GeantTrackGeo_v(4*maxPerBasket);
#ifndef VECCORE_CUDA
#ifdef USE_VECGEOM_NAVIGATOR
//  fRndm = &RNG::Instance();
  fRndm = new vecgeom::RNG; // what about the seed?
#elif USE_ROOT
  fRndm = new TRandom();
#endif
#endif
  fTransported = new GeantTrack_v(maxPerBasket, TrackDataMgr::GetInstance()->GetMaxDepth());
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTaskData::GeantTaskData(void *addr, size_t nthreads, int maxPerBasket, GeantPropagator *prop)
    : fPropagator(prop), fNthreads(nthreads)
{
  // Constructor
  char *buffer = (char*)addr;
  buffer += GeantTrack::round_up_align(sizeof(GeantTaskData));
  buffer = GeantTrack::round_up_align(buffer);

  fPath = VolumePath_t::MakeInstanceAt(TrackDataMgr::GetInstance()->GetMaxDepth(), (void*)buffer);
  fPathV = new VolumePath_t*[4*maxPerBasket];
  fNextpathV = new VolumePath_t*[4*maxPerBasket];
  fGeoTrack = GeantTrackGeo_v::MakeInstanceAt(buffer, 4*maxPerBasket);
  buffer += GeantTrackGeo_v::SizeOfInstance(4*maxPerBasket);
  buffer += VolumePath_t::SizeOfInstance(TrackDataMgr::GetInstance()->GetMaxDepth());
  buffer = GeantTrack::round_up_align(buffer);

  // Previous, the size was hard coded to 1024, '4' is a guess on the max number
  // of produced particles ...
  fTransported = GeantTrack_v::MakeInstanceAt(buffer, 4*maxPerBasket, TrackDataMgr::GetInstance()->GetMaxDepth());
  buffer += GeantTrack_v::SizeOfInstance(4*maxPerBasket, TrackDataMgr::GetInstance()->GetMaxDepth());

  fSizeInt = fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new (buffer) bool[fSizeBool];
  buffer += fSizeBool*sizeof(bool);
  fDblArray = new (buffer) double[fSizeDbl];
  buffer += fSizeDbl*sizeof(double);
  fIntArray = new (buffer) int[fSizeInt];
  buffer += fSizeInt*sizeof(int);

  fTrack = GeantTrack::MakeInstance();

#ifndef VECCORE_CUDA
#ifdef USE_VECGEOM_NAVIGATOR
  fRndm = &vecgeom::RNG::Instance();
#elif USE_ROOT
  fRndm = new TRandom();
#endif
#endif
//  fStat = new TrackStat(this);
}

//______________________________________________________________________________
GeantTaskData::~GeantTaskData()
{
// Destructor
  GeantTrack::ReleaseInstance(fTrack);
#ifndef VECCORE_CUDA
#ifndef USE_VECGEOM_NAVIGATOR
  delete fRndm;
#endif
#endif
  delete[] fBoolArray;
  delete[] fDblArray;
  delete[] fIntArray;
  delete [] fPathV;
  delete [] fNextpathV;
  delete fRndm;
  VolumePath_t::ReleaseInstance(fPath);
  delete fTransported;
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
void GeantTaskData::AttachPropagator(GeantPropagator *prop, int node)
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
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTaskData *GeantTaskData::MakeInstanceAt(void *addr, size_t nTracks, int maxPerBasket, GeantPropagator *prop)
{
   // GeantTrack MakeInstance based on a provided single buffer.
   return new (addr) GeantTaskData(addr, nTracks, maxPerBasket, prop);
}


//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t GeantTaskData::SizeOfInstance(size_t /*nthreads*/, int maxPerBasket)
{
   // @brief return the contiguous memory size needed to hold a GeantTrack_v

   const size_t bufSize = 5; // See constructor!

   size_t need = sizeof(GeantTaskData) // vecgeom::DevicePtr<Geant::cuda::GeantTaskData>::SizeOf()
      + GeantTrack::round_up_align(bufSize*maxPerBasket*(sizeof(bool)+sizeof(double)+sizeof(int)))
      + GeantTrack::round_up_align(VolumePath_t::SizeOfInstance(TrackDataMgr::GetInstance()->GetMaxDepth()))
      + GeantTrack_v::SizeOfInstance(4*maxPerBasket,TrackDataMgr::GetInstance()->GetMaxDepth());
   return GeantTrack::round_up_align(need);
}


#ifndef VECCORE_CUDA
//______________________________________________________________________________
GeantBasket *GeantTaskData::GetNextBasket()
{
  // Gets next free basket from the queue.
  if (fPool.empty())
    return nullptr;
  GeantBasket *basket = fPool.back();
  //  basket->Clear();
  fPool.pop_back();
  return basket;
}

//______________________________________________________________________________
void GeantTaskData::RecycleBasket(GeantBasket *b)
{
  // Recycle a basket.
  fPool.push_back(b);
}

//______________________________________________________________________________
void GeantTaskData::RecycleBasket(Basket *b)
{
  // Recycle a basket.
  fBPool.push_back(b);
}

//______________________________________________________________________________
int GeantTaskData::CleanBaskets(size_t ntoclean)
{
  // Clean a number of recycled baskets to free some memory
  GeantBasket *b;
  int ncleaned = 0;
  size_t ntodo = 0;
  if (ntoclean == 0)
    ntodo = fPool.size() / 2;
  else
    ntodo = Math::Min(ntodo, fPool.size());
  for (size_t i = 0; i < ntodo; i++) {
    b = fPool.back();
    delete b;
    ncleaned++;
    fPool.pop_back();
  }
  fToClean = false;
  //  Printf("Thread %d cleaned %d baskets", fTid, ncleaned);
  return ncleaned;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
GeantTrack &GeantTaskData::GetNewTrack()
{
  size_t index;
  if (fBlock->IsDistributed()) {
    fBlock = fPropagator->fTrackMgr->GetNewBlock();
    //printf("== New block: %d (%d) current=%d used=%d\n",
    //       fBlock->GetId(), fBlock->GetNode(), fBlock->GetCurrent(), fBlock->GetUsed());
    assert(fBlock->GetCurrent() == 0 && fBlock->GetUsed() == 0);
  }
  GeantTrack *track = fBlock->GetObject(index);
  track->Reset(*fTrack);
  track->SetBindex(index);
  return *track;
//  GeantTrack &track = fPropagator->fTrackMgr->GetTrack();
//  index = track.BIndex();
//  track.Reset(*fTrack);
//  track.SetBindex(index);
//  return track;
}

//______________________________________________________________________________
VECCORE_ATT_HOST_DEVICE
void GeantTaskData::ReleaseTrack(GeantTrack &track) {
  fPropagator->fTrackMgr->ReleaseTrack(track);
}

//______________________________________________________________________________
void GeantTaskData::InspectStages(int istage)
{
  Geant::Printf("** Thread %d: **", fTid);
  for (auto stage : fPropagator->fStages) {
    if (stage->GetId() == istage)
      Geant::Printf("*** -> %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());
    else
      Geant::Printf("***    %15s:  %d tracks", stage->GetName(), fStageBuffers[stage->GetId()]->size());    
  }
}


#endif

} // GEANT_IMPL_NAMESPACE
} // geant
