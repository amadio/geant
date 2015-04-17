#include "GeantThreadData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TArrayI.h"
#include "TGeoVolume.h"
#include "TRandom.h"

ClassImp(GeantThreadData)

//______________________________________________________________________________
GeantThreadData::GeantThreadData()
    : TObject(), fTid(-1), fNthreads(0), fMaxDepth(0), fSizeBool(0), fSizeDbl(0), 
      fVolume(0), fRndm(new TRandom()), fBoolArray(0), fDblArray(0), fTrack(0), 
      fPath(0), fBmgr(0), fPool() {
  // Constructor
  GeantPropagator *propagator = GeantPropagator::Instance();
  fNthreads = propagator->fNthreads;
  fMaxDepth = propagator->fMaxDepth;
  fSizeBool = fSizeDbl = 5 * propagator->fMaxPerBasket;
  fBoolArray = new Bool_t[fSizeBool];
  fDblArray = new Double_t[fSizeDbl];
  fPath = VolumePath_t::MakeInstance(fMaxDepth);
}

//______________________________________________________________________________
GeantThreadData::~GeantThreadData() {
  // Destructor
  //  delete fMatrix;
  delete fRndm;
  delete[] fBoolArray;
  delete[] fDblArray;
  VolumePath_t::ReleaseInstance(fPath);
}

//______________________________________________________________________________
Double_t *GeantThreadData::GetDblArray(Int_t size) {
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
Bool_t *GeantThreadData::GetBoolArray(Int_t size) {
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

//______________________________________________________________________________
GeantBasket *GeantThreadData::GetNextBasket() {
  // Gets next free basket from the queue.
  if (fPool.empty()) return 0;
  GeantBasket *basket = fPool.back();
  fPool.pop_back();
  return basket;
}  

//______________________________________________________________________________
void GeantThreadData::RecycleBasket(GeantBasket *b) {
  // Recycle a basket.
  fPool.push_back(b);
}

//______________________________________________________________________________
Int_t GeantThreadData::CleanBaskets(size_t ntoclean) {
  // Clean a number of recycled baskets to free some memory
  GeantBasket *b;
  Int_t ncleaned = 0;
  size_t ntodo = TMath::Min(fPool.size(), ntoclean);
  for (size_t i=0; i<ntodo; i++) {
    b = fPool.back();
    delete b;
    ncleaned++;
    fPool.pop_back();
  }
  return ncleaned;
}
  
