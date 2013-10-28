#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TMutex.h"
//==============================================================================
// GeantScheduler - dispatcher running in a single thread. Collects tracks
// from all threads via an input queue and fills baskets corresponding to each
// volume, which are then injected in the main work queue.
//==============================================================================

class concurrent_queue;
class GeantTrack_v;

//______________________________________________________________________________
class GeantScheduler : public TObject {
protected:
   Int_t                fNvolumes;            // Number of active volumes in the geometry
   Int_t                fNpriority;           // Number of priority baskets held
   GeantBasketMgr     **fBasketMgr;           // Array of basket managers
   Int_t                fPriorityRange[2];    // Prioritized events
   
public:
   GeantScheduler();
   virtual ~GeantScheduler();
   Int_t                AddTracks(GeantTrack_v &tracks);
   void                 CreateBaskets();
   Int_t                CollectPrioritizedTracks();
   Int_t                GetNpriority() const {return fNpriority;}
   void                 SetPriorityRange(Int_t min, Int_t max) {fPriorityRange[0]=min; fPriorityRange[1]=max;}
   Int_t                FlushPriorityBaskets();
   Int_t                GarbageCollect();
   
   ClassDef(GeantScheduler, 1)      // Main basket scheduler
};   
#endif
