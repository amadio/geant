#ifndef EVENT_HIT
#define EVENT_HIT

#include "Rtypes.h"

#include "tbb/atomic.h"

//______________________________________________________________________________
struct GeantEvent {
   Int_t  event;        // event number
   Int_t  slot;         // fixed slot number
   tbb::atomic<Int_t>  ntracks;      // number of tracks
   tbb::atomic<Int_t>  ndone;        // number of done tracks

   GeantEvent() : event(0), slot(0), ntracks(), ndone() {}
   GeantEvent(Int_t ntr) : event(0), slot(0), ntracks(), ndone() { ntracks=ntr; }
   ~GeantEvent() {}
   
   void               AddTrack();
   void               Reset() {ntracks = ndone = 0;}
   void               StopTrack();
   Bool_t             Transported() const {return ((ntracks>0) && (ntracks==ndone));}
   Int_t              TracksLeft() const {return ntracks-ndone;}
};

//______________________________________________________________________________
struct GeantHit {
   Double_t x;        // X position
   Double_t y;        // Y position
   Double_t z;        // Z position
   Double_t eloss;    // Energy loss
   Int_t    event;    // Event number
   
   GeantHit() : x(0.), y(0.), z(0.), eloss(0.), event(-1) {}
   GeantHit(Double_t xh, Double_t yh, Double_t zh, Double_t elossh, Int_t eventh)
              : x(xh), y(yh), z(zh), eloss(elossh), event(eventh) {}
            
};

#endif // EVENT_HIT
