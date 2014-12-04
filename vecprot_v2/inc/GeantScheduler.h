#ifndef GEANT_SCHEDULER
#define GEANT_SCHEDULER

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifdef __STAT_DEBUG
#include "GeantTrackStat.h"
#endif   

#if __cplusplus >= 201103L
#include <atomic>
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
   GeantBasketMgr      *fGarbageCollector;    // Garbage collector manager
#if __cplusplus >= 201103L
   std::atomic_int     *fNtracks;             //[fNvolumes] Number of tracks per volume
#else
   Int_t               *fNtracks;             //[fNvolumes] Number of tracks per volume
#endif
   Int_t                fPriorityRange[2];    // Prioritized events
#ifdef __STAT_DEBUG
   GeantTrackStat       fPStat;  //! Statistics for the pending tracks
   GeantTrackStat       fQStat;  //! Statistics for the queued tracks
   GeantTrackStat       fTStat;  //! Statistics for the transported tracks
#endif   

private:
   GeantScheduler(const GeantScheduler&);
   GeantScheduler& operator=(const GeantScheduler&);
  
public:
   GeantScheduler();
   virtual ~GeantScheduler();
   Int_t                AddTrack(GeantTrack &track);
   Int_t                AddTracks(GeantBasket *output, Int_t &ntot, Int_t &nnew, Int_t &nkilled);
   void                 AdjustBasketSize();
   void                 CreateBaskets();
   Int_t                CollectPrioritizedTracks();
   GeantBasketMgr     **GetBasketManagers() const {return fBasketMgr;}
   GeantBasketMgr      *GetGarbageCollector() const {return fGarbageCollector;}
#if __cplusplus >= 201103L
   Int_t                GetNtracks(Int_t ib) {return fNtracks[ib].load();}
#else   
   Int_t                GetNtracks(Int_t ib) {return fNtracks[ib];}
#endif
   Int_t                GetNpriority() const {return fNpriority;}
   Int_t                GetNvolumes() const  {return fNvolumes;}
   void                 SetPriorityRange(Int_t min, Int_t max) {fPriorityRange[0]=min; fPriorityRange[1]=max;}
#ifdef __STAT_DEBUG
   GeantTrackStat      &GetPendingStat() {return fPStat;}
   GeantTrackStat      &GetQueuedStat()  {return fQStat;}
   GeantTrackStat      &GetTransportStat()  {return fTStat;}
#endif   
   Int_t                FlushPriorityBaskets();
   Int_t                GarbageCollect();
   void                 PrintSize() const;
   size_t               Sizeof() const;
      
   ClassDef(GeantScheduler, 1)      // Main basket scheduler
};   
#endif
