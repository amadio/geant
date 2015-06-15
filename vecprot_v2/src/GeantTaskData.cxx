#include "GeantTaskData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TArrayI.h"
#include "TGeoVolume.h"
#include "TRandom.h"

ClassImp(GeantTaskData)

//______________________________________________________________________________
GEANT_CUDA_DEVICE_CODE
GeantTaskData::GeantTaskData(Int_t nthreads, Int_t maxDepth, Int_t maxPerBasket)
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(0), fDblArray(0), fTrack(0,maxDepth),
      fPath(0), fBmgr(0), fPool() {
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

//______________________________________________________________________________
GeantTaskData::GeantTaskData()
    : fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), fToClean(false),
      fVolume(0), fRndm(nullptr), fBoolArray(0), fDblArray(0), fTrack(0),
      fPath(0), fBmgr(0), fPool() {
  // Constructor
  GeantPropagator *propagator = GeantPropagator::Instance();
  fNthreads = propagator->fNthreads;
  fMaxDepth = propagator->fMaxDepth;
  fSizeBool = fSizeDbl = 5 * propagator->fMaxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
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
  VolumePath_t::ReleaseInstance(fPath);
}

//______________________________________________________________________________
Double_t *GeantTaskData::GetDblArray(Int_t size) {
  // Return memory storage for an array of doubles of at least "size" length which
  // is thread-specific
  if (size < fSizeDbl)
    return fDblArray;
  Double_t *array = new Double_t[size];
  delete[] fDblArray;
  fDblArray = array;
  fSizeDbl = size;
  return fDblArray;
}

//______________________________________________________________________________
Bool_t *GeantTaskData::GetBoolArray(Int_t size) {
  // Return memory storage for an array of bool of at least "size" length which
  // is thread-specific
  if (size < fSizeBool)
    return fBoolArray;
  Bool_t *array = new Bool_t[size];
  delete[] fBoolArray;
  fBoolArray = array;
  fSizeBool = size;
  return fBoolArray;
}

#ifndef GEANT_NVCC
//______________________________________________________________________________
GeantBasket *GeantTaskData::GetNextBasket() {
  // Gets next free basket from the queue.
  if (fPool.empty()) return 0;
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
