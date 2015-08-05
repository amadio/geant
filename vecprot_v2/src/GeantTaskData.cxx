#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"
#include "Geant/Typedefs.h"

#ifdef USE_ROOT
#include "TRandom.h"
#endif

#include "base/SOA3D.h"

using std::min;

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::GeantTaskData(int nthreads, int maxDepth, int maxPerBasket)
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false), fVolume(0), fRndm(nullptr),
      fBoolArray(nullptr), fDblArray(nullptr), fTrack(0, maxDepth), fPath(0), fBmgr(0), fPool(),
      fSOA3Dworkspace1(new vecgeom::SOA3D<vecgeom::Precision>(5 * maxPerBasket)),
      fSOA3Dworkspace2(new vecgeom::SOA3D<vecgeom::Precision>(5 * maxPerBasket)), fSizeInt(5 * maxPerBasket),
      fIntArray(new int[fSizeInt]) {
  // Constructor
  fNthreads = nthreads;
  fMaxDepth = maxDepth;
  fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new bool[fSizeBool];
  fDblArray = new double[fSizeDbl];
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
#ifndef GEANT_NVCC
#ifdef USE_ROOT
  fRndm = new TRandom();
#else
  fRndm = &RNG::Instance();
#endif
#endif
}

//______________________________________________________________________________
GeantTaskData::GeantTaskData()
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false), fVolume(0), fRndm(nullptr),
      fBoolArray(nullptr), fDblArray(nullptr), fTrack(0), fPath(0), fBmgr(0), fPool(), fSOA3Dworkspace1(),
      fSOA3Dworkspace2(), fSizeInt(0), fIntArray(nullptr) {
  // Constructor
  GeantPropagator *propagator = GeantPropagator::Instance();
  fNthreads = propagator->fNthreads;
  fMaxDepth = propagator->fMaxDepth;
  fSizeBool = fSizeDbl = fSizeInt = 5 * propagator->fMaxPerBasket;
  fBoolArray = new bool[fSizeBool];
  fDblArray = new double[fSizeDbl];
  fIntArray = new int[fSizeInt];
  fSOA3Dworkspace1 = new vecgeom::SOA3D<double>(fSizeInt);
  fSOA3Dworkspace2 = new vecgeom::SOA3D<double>(fSizeInt);
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
#ifdef USE_ROOT
  fRndm = new TRandom();
#else
  fRndm = &RNG::Instance();
#endif
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::~GeantTaskData() {
// Destructor
//  delete fMatrix;
#ifndef GEANT_NVCC
  delete fRndm;
#endif
  delete[] fBoolArray;
  delete[] fDblArray;
  delete[] fIntArray;
  delete fSOA3Dworkspace1;
  delete fSOA3Dworkspace2;
  VolumePath_t::ReleaseInstance(fPath);
}

#ifndef GEANT_NVCC
//______________________________________________________________________________
GeantBasket *GeantTaskData::GetNextBasket() {
  // Gets next free basket from the queue.
  if (fPool.empty())
    return nullptr;
  GeantBasket *basket = fPool.back();
  //  basket->Clear();
  fPool.pop_back();
  return basket;
}

//______________________________________________________________________________
void GeantTaskData::RecycleBasket(GeantBasket *b) {
  // Recycle a basket.
  fPool.push_back(b);
}

//______________________________________________________________________________
int GeantTaskData::CleanBaskets(size_t ntoclean) {
  // Clean a number of recycled baskets to free some memory
  GeantBasket *b;
  int ncleaned = 0;
  size_t ntodo = 0;
  if (ntoclean == 0)
    ntodo = fPool.size() / 2;
  else
    ntodo = min<int>(ntodo, fPool.size());
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

#endif

} // GEANT_IMPL_NAMESPACE
} // geant
