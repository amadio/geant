#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifdef __STAT_DEBUG
#include "GeantTrackStat.h"
#endif   

#include "TMutex.h"
//==============================================================================
// GeantScheduler - dispatcher running in a single thread. Collects tracks
// from all threads via an input queue and fills baskets corresponding to each
// volume, which are then injected in the main work queue.
//==============================================================================

class concurrent_queue;
class GeantTrack;
class GeantBasket;
class GeantBasketMgr;

//______________________________________________________________________________
class GeantScheduler : public TObject {
protected:
   Int_t                fNvolumes;            // Number of active volumes in the geometry
   Int_t                fNpriority;           // Number of priority baskets held
   GeantBasketMgr     **fBasketMgr;           // Array of basket managers
   Int_t               *fNtracks;             //[fNvolume] Number of tracks per volume
   Int_t                fPriorityRange[2];    // Prioritized events
#ifdef __STAT_DEBUG
   GeantTrackStat       fPStat;  //! Statistics for the pending tracks
   GeantTrackStat       fQStat;  //! Statistics for the queued tracks
   GeantTrackStat       fTStat;  //! Statistics for the transported tracks
#endif   
   
public:
   GeantScheduler();
   virtual ~GeantScheduler();
   Int_t                AddTrack(const GeantTrack &track);
   Int_t                AddTracks(GeantBasket *output, Int_t &ntot, Int_t &nnew, Int_t &nkilled);
   void                 AdjustBasketSize();
   void                 CreateBaskets();
   Int_t                CollectPrioritizedTracks();
   Int_t                GetNpriority() const {return fNpriority;}
   void                 SetPriorityRange(Int_t min, Int_t max) {fPriorityRange[0]=min; fPriorityRange[1]=max;}
#ifdef __STAT_DEBUG
   GeantTrackStat      &GetPendingStat() {return fPStat;}
   GeantTrackStat      &GetQueuedStat()  {return fQStat;}
   GeantTrackStat      &GetTransportStat()  {return fTStat;}
#endif   
   Int_t                FlushPriorityBaskets();
   Int_t                GarbageCollect();
   
   ClassDef(GeantScheduler, 1)      // Main basket scheduler
};   
#endif
