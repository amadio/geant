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
   Int_t                fNvolumes;            // Number of active volumes in the geometry
   Int_t                fNpriority;           // Number of priority baskets held
   GeantBasket        **fBaskets;             // Array of baskets to be filled
   GeantBasket        **fPriorityBaskets;     // Array of priority baskets

   TMutex the_lock;

public:
   GeantMainScheduler();
   GeantMainScheduler(Int_t nvolumes);
   virtual ~GeantMainScheduler();

   static GeantMainScheduler* fgInstance;
   static GeantMainScheduler* Instance(Int_t nvolumes=0);

   Int_t                AddTrack(Int_t itrack, Int_t ibasket, Bool_t* pushedPriority);
   Int_t                GetNpriority() const {return fNpriority;}

   Int_t                FlushNormalBaskets(Int_t threshold=0);
   Int_t                FlushPriorityBaskets();

   ClassDef(GeantMainScheduler, 1)      // Main basket scheduler
};

#endif // GEANT_MAINSCHEDULER
