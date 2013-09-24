#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TMutex.h"
//==============================================================================
// GeantMainScheduler - dispatcher running in a single thread. Collects tracks
// from all threads via an input queue and fills baskets corresponding to each
// volume, which are then injected in the main work queue.
//==============================================================================

class concurrent_queue;
class GeantTrack_v;

//______________________________________________________________________________
class GeantMainScheduler : public TObject {
protected:
   Int_t                fNvolumes;            // Number of active volumes in the geometry
   Int_t                fNpriority;           // Number of priority baskets held
   GeantBasket        **fBaskets;             // Array of baskets to be filled
   GeantBasket        **fPriorityBaskets;     // Array of priority baskets
   Int_t                fPriorityRange[2];    // Prioritized events
   concurrent_queue    *feeder_queue;         // Main feeder
   concurrent_queue    *empty_queue;          // Queue with empty baskets
   concurrent_queue    *collector_queue;      // Queue for collecting tracks
   
public:
   GeantMainScheduler();
   GeantMainScheduler(Int_t nvolumes);
   virtual ~GeantMainScheduler();
   Int_t                AddTracks(GeantTrack_v &tracks);
   Int_t                AddTrack(Int_t itrack, Int_t ibasket);
   Int_t                GetNpriority() const {return fNpriority;}
   void                 SetPriorityRange(Int_t min, Int_t max) {fPriorityRange[0]=min; fPriorityRange[1]=max;}
   Int_t                FlushBaskets(Int_t threshold=0);   
   Int_t                FlushPriorityBaskets();
   
   ClassDef(GeantMainScheduler, 1)      // Main basket scheduler
};   
#endif
