#ifndef GEANT_MAINSCHEDULER
#define GEANT_MAINSCHEDULER

//==============================================================================
// GeantMainScheduler - dispatcher running in a single thread. Collects tracks
// from all threads via an input queue and fills baskets corresponding to each
// volume, which are then injected in the main work queue.
//==============================================================================

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

class GeantBasket;

//______________________________________________________________________________
class GeantMainScheduler : public TObject {
protected:
   int                fNvolumes;            // Number of active volumes in the geometry
   GeantBasket        **fBaskets;             // Array of baskets to be filled
   GeantBasket        **fPriorityBaskets;     // Array of priority baskets

   TMutex the_lock;

public:
   GeantMainScheduler();
   GeantMainScheduler(int nvolumes);
   virtual ~GeantMainScheduler();

   static GeantMainScheduler* fgInstance;
   static GeantMainScheduler* Instance(int nvolumes=0);

   int                AddTrack(int itrack, int ibasket, bool* pushedPriority);

   int                FlushNormalBaskets();
   int                FlushPriorityBaskets();

   ClassDef(GeantMainScheduler, 1)      // Main basket scheduler
};

#endif // GEANT_MAINSCHEDULER
