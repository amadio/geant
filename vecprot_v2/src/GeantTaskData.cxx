#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "Basket.h"
#include "StackLikeBuffer.h"
#include "GeantPropagator.h"
#include "TrackManager.h"
#include "GeantTrackGeo.h"
#include "Geant/Typedefs.h"
#include "Geant/Error.h"
#include "SimulationStage.h"
#include "TrackStat.h"

#ifdef USE_ROOT
#include "TRandom.h"
#endif

using std::min;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GeantTaskData::GeantTaskData(size_t nthreads, int maxDepth, int maxPerBasket)
    : fPropagator(nullptr), fTid(-1), fNode(0), fNthreads(nthreads), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(nullptr), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0, maxDepth),
      fPath(nullptr), fBmgr(nullptr), fReused(nullptr), fImported(nullptr), fStackBuffer(nullptr), fPool(),
      fSizeInt(5 * maxPerBasket), fIntArray(new int[fSizeInt]), fTransported(nullptr), fTransported1(maxPerBasket), fNkeepvol(0),
      fNsteps(0), fNsnext(0), fNphys(0), fNmag(0), fNpart(0), fNsmall(0), fNcross(0), fPhysicsData(nullptr)
{
  // Constructor
  fNthreads = nthreads;
  fMaxDepth = maxDepth;
  fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new bool[fSizeBool];
  fDblArray = new double[fSizeDbl];
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
  fPathV = new VolumePath_t*[maxPerBasket];
  fNextpathV = new VolumePath_t*[maxPerBasket];
  fGeoTrack = new GeantTrackGeo_v(maxPerBasket);
#ifndef VECCORE_CUDA
#ifdef USE_VECGEOM_NAVIGATOR
//  fRndm = &RNG::Instance();
  fRndm = new vecgeom::RNG; // what about the seed?
#elif USE_ROOT
  fRndm = new TRandom();
#endif
#endif
  fTransported = new GeantTrack_v(maxPerBasket, maxDepth);
  fShuttleBasket = new Basket(1000, 0, -1);
  fBvector = new Basket(256, 0, -1);
  fStat = new TrackStat(this);
  for (int i=0; i<=int(kSteppingActionsStage); ++i)
    fStageBuffers.push_back(new Basket(1000, 0, -1));
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTaskData::GeantTaskData(void *addr, size_t nthreads, int maxDepth, int maxPerBasket, GeantPropagator *prop /* = nullptr */)
    : fPropagator(prop), fTid(-1), fNode(0), fNthreads(nthreads), fMaxDepth(maxDepth), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(nullptr), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0, maxDepth),
      fPath(nullptr), fBmgr(nullptr), fReused(nullptr), fImported(nullptr), fStackBuffer(nullptr), fPool(),
      fSizeInt( 5*maxPerBasket ), fIntArray( nullptr ), fTransported(nullptr), fNkeepvol(0),
      fNsteps(0), fNsnext(0), fNphys(0), fNmag(0), fNpart(0), fNsmall(0), fNcross(0), fPhysicsData(nullptr)
{
  // Constructor
  char *buffer = (char*)addr;
  buffer += GeantTrack::round_up_align(sizeof(GeantTaskData));
  const size_t nElements = 5; // See other constructor!

  buffer = GeantTrack::round_up_align(buffer);

  fPath = VolumePath_t::MakeInstanceAt(fMaxDepth,(void*)buffer);
  fPathV = new VolumePath_t*[maxPerBasket];
  fNextpathV = new VolumePath_t*[maxPerBasket];
  fGeoTrack = GeantTrackGeo_v::MakeInstanceAt(buffer, 4*maxPerBasket);
  buffer += GeantTrackGeo_v::SizeOfInstance(4*maxPerBasket);
  buffer += VolumePath_t::SizeOfInstance(fMaxDepth);
  buffer = GeantTrack::round_up_align(buffer);

  // Previous, the size was hard coded to 1024, '4' is a guess on the max number
  // of produced particles ...
  fTransported = GeantTrack_v::MakeInstanceAt(buffer, 4*maxPerBasket, fMaxDepth);
  buffer += GeantTrack_v::SizeOfInstance(4*maxPerBasket, fMaxDepth);

  fSizeInt = fSizeBool = fSizeDbl = nElements * maxPerBasket;
  fBoolArray = new (buffer) bool[fSizeBool];
  buffer += fSizeBool*sizeof(bool);
  fDblArray = new (buffer) double[fSizeDbl];
  buffer += fSizeDbl*sizeof(double);
  fIntArray = new (buffer) int[fSizeInt];
  buffer += fSizeInt*sizeof(int);


#ifndef VECCORE_CUDA
#ifdef USE_VECGEOM_NAVIGATOR
  fRndm = &vecgeom::RNG::Instance();
#elif USE_ROOT
  fRndm = new TRandom();
#endif
#endif
  fShuttleBasket = new Basket(1000, 0, -1);
  fBvector = new Basket(256, 0, -1);
  fStat = new TrackStat(this);
  for (int i=0; i<=int(kSteppingActionsStage); ++i)
    fStageBuffers.push_back(new Basket(1000, 0, prop->fNuma));
}

//______________________________________________________________________________
GeantTaskData::~GeantTaskData()
{
// Destructor
//  delete fMatrix;
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
  for (int i=0; i<=int(kSteppingActionsStage); ++i)
    delete fStageBuffers[i];
}

//______________________________________________________________________________
VECCORE_ATT_DEVICE
GeantTaskData *GeantTaskData::MakeInstanceAt(void *addr, size_t nTracks, int maxdepth, int maxPerBasket, GeantPropagator *prop)
{
   // GeantTrack MakeInstance based on a provided single buffer.
   return new (addr) GeantTaskData(addr, nTracks, maxdepth, maxPerBasket, prop);
}


//______________________________________________________________________________
VECCORE_ATT_DEVICE
size_t GeantTaskData::SizeOfInstance(size_t /*nthreads*/, int maxDepth, int maxPerBasket)
{
   // @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth

   const size_t bufSize = 5; // See constructor!

   size_t need = sizeof(GeantTaskData) // vecgeom::DevicePtr<Geant::cuda::GeantTaskData>::SizeOf()
      + GeantTrack::round_up_align(bufSize*maxPerBasket*(sizeof(bool)+sizeof(double)+sizeof(int)))
      + GeantTrack::round_up_align(VolumePath_t::SizeOfInstance(maxDepth))
      + GeantTrack_v::SizeOfInstance(4*maxPerBasket,maxDepth);
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
Basket *GeantTaskData::GetFreeBasket()
{
  // Gets next free basket from the queue.
  if (fBPool.empty())
    return ( new Basket(fPropagator->fConfig->fMaxPerBasket) );
  Basket *basket = fBPool.back();
  //  basket->Clear();
  fBPool.pop_back();
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
  track->Clear();
  track->fBindex = index;
  track->fMaxDepth = fMaxDepth;
  return *track;
  
//  return ( fPropagator->fTrackMgr->GetTrack() );
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
