#ifndef GEANT_PARTICLEBUFFER
#define GEANT_PARTICLEBUFFER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantVolumeBasket;

// Buffer of pending particles. Independent storage per thread. Flushing should
// be done by a single thread.

//______________________________________________________________________________
class GeantParticleBuffer : public TObject {
protected:
   Int_t               fNthreads;              // number of threads
   Int_t               fCrtFlushed;            // current thread being flushed
   Int_t               fMaxPending;            // max number of pending tracks per thread
   Int_t              *fNpending;              //[fNthreads] Number of pending particles per thread
   Int_t             **fTrackInd;              //[fNthreads*fMaxPending]
   GeantVolumeBasket ***fBaskets;               //[fNthreads*fMaxPending]

public:
   GeantParticleBuffer() : TObject(), fNthreads(0), fCrtFlushed(0), fMaxPending(0), fNpending(0), fTrackInd(0), fBaskets(0) {}
   GeantParticleBuffer(Int_t nthreads, Int_t max_pending);
   ~GeantParticleBuffer();
   
   void                AddPendingTrack(Int_t itrack, GeantVolumeBasket *basket, Int_t tid);
   Bool_t              FlushBaskets();
   void                Reset() {fCrtFlushed = 0;}
   Bool_t              IsEmpty() const;
   
   ClassDef(GeantParticleBuffer, 0)
};   
#endif
