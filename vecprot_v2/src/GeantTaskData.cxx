#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TArrayI.h"
#include "TGeoVolume.h"
#include "TRandom.h"

#include "base/SOA3D.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::GeantTaskData(size_t nthreads, Int_t maxDepth, Int_t maxPerBasket)
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0,maxDepth),
      fPath(0), fBmgr(0), fPool(),
      fSOA3Dworkspace1( new vecgeom::SOA3D<vecgeom::Precision>(5*maxPerBasket) ),
      fSOA3Dworkspace2( new vecgeom::SOA3D<vecgeom::Precision>(5*maxPerBasket) ),
      fSizeInt( 5*maxPerBasket ),
      fIntArray( new int[fSizeInt] )
      {
  // Constructor
  fNthreads = nthreads;
  fMaxDepth = maxDepth;
  fSizeBool = fSizeDbl = 5 * maxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
#ifndef GEANT_NVCC
  fRndm = new TRandom();
#endif
}

template <typename T>
GEANT_CUDA_DEVICE_CODE
size_t soaSizeOfInstance(size_t nElements)
{
   return GeantTrack_v::round_up_align(sizeof(vecgeom::SOA3D<T>)+3*sizeof(T)*nElements);
}

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::GeantTaskData(void *addr, size_t nthreads, Int_t maxDepth, Int_t maxPerBasket)
    : fTid(-1), fNthreads(nthreads), fMaxDepth(maxDepth), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(nullptr), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0,maxDepth),
      fPath(nullptr), fBmgr(nullptr), fPool(),
      fSOA3Dworkspace1(nullptr),
      fSOA3Dworkspace2(nullptr),
      fSizeInt( 5*maxPerBasket ),
      fIntArray( nullptr )
      {
  // Constructor
  char *buffer = (char*)addr;
  buffer += sizeof(GeantTaskData);
  const size_t nElements = 5; // See other constructor!

  vecgeom::Precision *x = (vecgeom::Precision*)(buffer + sizeof(vecgeom::Precision));
  vecgeom::Precision *y = x+(nElements*maxPerBasket);
  vecgeom::Precision *z = y+(nElements*maxPerBasket);
  fSOA3Dworkspace1 = new (buffer) vecgeom::SOA3D<vecgeom::Precision>(x,y,z,nElements*maxPerBasket);
  buffer += soaSizeOfInstance<vecgeom::Precision>(nElements*maxPerBasket);

  x = (vecgeom::Precision*)(buffer + sizeof(vecgeom::Precision));
  y = x+(nElements*maxPerBasket);
  z = y+(nElements*maxPerBasket);
  fSOA3Dworkspace2 = new (buffer) vecgeom::SOA3D<vecgeom::Precision>(x,y,z,nElements*maxPerBasket);
  buffer += soaSizeOfInstance<vecgeom::Precision>(nElements*maxPerBasket);

  fPath = VolumePath_t::MakeInstanceAt(fMaxDepth,(void*)buffer);
  buffer += VolumePath_t::SizeOfInstance(fMaxDepth);

  fSizeInt = fSizeBool = fSizeDbl = nElements * maxPerBasket;
  fBoolArray = new (buffer) Bool_t[fSizeBool];
  buffer += fSizeBool*sizeof(bool);
  fDblArray = new (buffer) Double_t[fSizeDbl];
  buffer += fSizeDbl*sizeof(double);
  fIntArray = new (buffer) int[fSizeInt];
  buffer += fSizeInt*sizeof(int);

#ifndef GEANT_NVCC
  fRndm = new TRandom();
#endif
}

//______________________________________________________________________________
GeantTaskData::GeantTaskData()
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(nullptr), fDblArray(nullptr), fTrack(0),
      fPath(0), fBmgr(0), fPool(),
      fSOA3Dworkspace1(),
      fSOA3Dworkspace2(),
      fSizeInt(0),
      fIntArray(nullptr) {
  // Constructor
  GeantPropagator *propagator = GeantPropagator::Instance();
  fNthreads = propagator->fNthreads;
  fMaxDepth = propagator->fMaxDepth;
  fSizeBool = fSizeDbl = fSizeInt = 5 * propagator->fMaxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
  fIntArray = new int[fSizeInt];
  fSOA3Dworkspace1 = new vecgeom::SOA3D<double>( fSizeInt );
  fSOA3Dworkspace2 = new vecgeom::SOA3D<double>( fSizeInt );
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
  fRndm = new TRandom();
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

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData *GeantTaskData::MakeInstanceAt(void *addr, size_t nTracks, Int_t maxdepth, Int_t maxPerBasket)
{
   // GeantTrack MakeInstance based on a provided single buffer.
   return new (addr) GeantTaskData(addr, nTracks, maxdepth, maxPerBasket);
}


//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
size_t GeantTaskData::SizeOfInstance(size_t nthreads, Int_t maxDepth, Int_t maxPerBasket)
{
   // @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth

   const size_t bufSize = 5; // See constructor!
   const size_t soaSize = soaSizeOfInstance<vecgeom::Precision>(bufSize*maxPerBasket);

   size_t need = sizeof(GeantTaskData) // vecgeom::DevicePtr<Geant::cuda::GeantTaskData>::SizeOf()
      + bufSize*maxPerBasket*(sizeof(bool)+sizeof(double)+sizeof(int))
      + VolumePath_t::SizeOfInstance(maxDepth)
      + 2*soaSize;
   return GeantTrack_v::round_up_align(need);
}


#ifndef GEANT_NVCC
GeantBasket *GeantTaskData::GetNextBasket() {
  // Gets next free basket from the queue.
  if (fPool.empty()) return nullptr;
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
Int_t GeantTaskData::CleanBaskets(size_t ntoclean) {
  // Clean a number of recycled baskets to free some memory
  GeantBasket *b;
  Int_t ncleaned = 0;
  size_t ntodo = 0;
  if (ntoclean==0) 
    ntodo = fPool.size()/2;
  else 
    ntodo = TMath::Min(fPool.size(), ntoclean);
  for (size_t i=0; i<ntodo; i++) {
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

