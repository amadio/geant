#ifndef EVENT_HIT
#define EVENT_HIT

#include "Rtypes.h"

#include "tbb/atomic.h"

//______________________________________________________________________________
struct GeantEvent {
   int  event;        // event number
   int  slot;         // fixed slot number
   tbb::atomic<int>  ntracks;      // number of tracks
   tbb::atomic<int>  ndone;        // number of done tracks

   GeantEvent() : event(0), slot(0), ntracks(), ndone() {}
   GeantEvent(int ntr) : event(0), slot(0), ntracks(), ndone() { ntracks=ntr; }
   ~GeantEvent() {}
   
   void               AddTrack();
   void               Reset() {ntracks = ndone = 0;}
   void               StopTrack();
   Bool_t             Transported() const {return ((ntracks>0) && (ntracks==ndone));}
   int              TracksLeft() const {return ntracks-ndone;}
};

//______________________________________________________________________________
struct GeantHit {
   double x;        // X position
   double y;        // Y position
   double z;        // Z position
   double eloss;    // Energy loss
   int    event;    // Event number
   
   GeantHit() : x(0.), y(0.), z(0.), eloss(0.), event(-1) {}
   GeantHit(double xh, double yh, double zh, double elossh, int eventh)
              : x(xh), y(yh), z(zh), eloss(elossh), event(eventh) {}
            
};

#endif // EVENT_HIT
