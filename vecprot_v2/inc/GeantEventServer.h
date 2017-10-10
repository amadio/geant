#ifndef GEANT_EVENT_SERVER_H
#define GEANT_EVENT_SERVER_H

#include <atomic>

#include <list>

#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantConfig.h"
#include "priority_queue.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PrimaryGenerator;
class GeantTaskData;
class GeantEvent;
class Basket;
class StackLikeBuffer;

//------------------------------------------------------------------------------
//
// A concurrent event server supporting:
// - Non-concurrent filling of tracks from input events
// - Concurrent retrieval of tracks from arbitrary number of clients
//
//  |---ev.0-----|----ev.1---|---ev.2---|--->    fTracks
//  ||||||||||||||||---------|||                 Tracks actually copied
//  fNevents                               ->    Number of events buffered
//  fStartEvt[0] fStartEvt[1] fStartEvt[2] ->    Start track index for each evt.
//  fNtracks[0]  fNtracks[1]  fNtracks[2]  ->    Number of tracks per event
//  fNFilled[0]  fNFilled[1]  fNfilled[2]  ->    Number of filled tracks per evt.
//  fNDispatched[0] ...                    ->    Number of dispatched per evt.
//------------------------------------------------------------------------------

class GeantEventServer
{
public:

enum ErrorType {
  kNoerr   = 0x00,    // normal operation
  kCSlots  = 0x01,    // contention on slots
  kCEvents = 0x02,    // contention on queued events
  kDone    = 0x04     // all events dispatched
};

using queue_slots = mpmc_bounded_queue<size_t>;
using queue_events = mpmc_bounded_queue<GeantEvent*>;

private:
  int fNevents = 0;                    /** Number of events to be filled */
  int fNprimaries = 0;                 /** Number of primaries served */
  int fNactiveMax = 0;                 /** Maximum number of active events (buffer size)*/
  int fNtracksInit = 0;                /** Initial number of tracks in the buffer. */
  int fNbasketsInit = 0;               /** Initial number of baskets to be served */
  std::atomic_int fNactive;            /** Number of deployed events */
  std::atomic_int fNserved;            /** Number of baskets served */
  std::atomic_int fLastActive;         /** Last activated event */
  std::atomic_int fCurrentEvent;       /** Current event being served */
  std::atomic_int fNload;              /** Last load event in the server */
  std::atomic_int fNstored;            /** Number of stored events in the server */
  std::atomic_int fNcompleted;         /** Number of completed events */
  std::atomic_flag fGenLock;           /** Generator lock */
  std::atomic<GeantEvent*> fEvent;     /** Current event being distributed */
  GeantRunManager *fRunMgr = nullptr;  /** Run manager */
  bool fEventsServed = false;          /** All events served */
  bool fDone = false;                  /** All events transported */
  bool fHasTracks = false;             /** Server has tracks to dispatch */
  bool fInitialPhase = false;          /** Server in initial dispatch phase */
  GeantEvent** fEvents = nullptr;      /** Events to be dispatched */
  int  fBindex = 0;                    /** Basket manager index */
  queue_slots fFreeSlots;              /** Queue of free event slots */
  queue_events fPendingEvents;         /** Queue of pending events */
  queue_events fDoneEvents;            /** Queue of transported events */

protected:
  GeantTrack *GetNextTrack(unsigned int &error);

public:
  GeantEventServer(int nactive_max, GeantRunManager *runmgr);
  ~GeantEventServer();

  GEANT_FORCE_INLINE
  size_t  AdjustSize(size_t size) const
  { size_t i = 1, new_size = 2; while ((size >> i++) > 0) new_size *= 2; return new_size; }

// Accessors
  GEANT_FORCE_INLINE
  int  GetNevents() const { return fNevents; }

  GEANT_FORCE_INLINE
  int  GetNprimaries() const { return fNprimaries; }

  GEANT_FORCE_INLINE
  int  GetNbasketsInit() const { return fNbasketsInit; }

  GEANT_FORCE_INLINE
  int  GetNactiveMax() const { return fNactiveMax; }

  GEANT_FORCE_INLINE
  int  GetNactive() const { return fNactive.load(); }

  GEANT_FORCE_INLINE
  int  GetCurrent() const { return fCurrentEvent.load(); }

  GEANT_FORCE_INLINE
  int  GetNstored() const { return fNstored.load(); }
  
  GEANT_FORCE_INLINE
  GeantEvent *GetEvent(int slot) { return fEvents[slot]; }

  GEANT_FORCE_INLINE
  GeantEvent *FindEvent(int event) {
    // trying to avoid concurrent map, but this could be smarter
    for (int i=0; i<fNactiveMax; ++i) {
      if (fEvents[i]->GetEvent() == event) return fEvents[i];
    }
    return nullptr;
  }

  GEANT_FORCE_INLINE
  int GetBindex() { return fBindex; }

  GEANT_FORCE_INLINE
  bool EventsServed() const { return fEventsServed; }

  GEANT_FORCE_INLINE
  bool HasTracks() const { return fHasTracks; }

  GEANT_FORCE_INLINE
  bool IsInitialPhase() const { return fInitialPhase; }

  int FillBasket(GeantTrack_v &tracks, int ntracks, unsigned int &error);

  int FillBasket(Basket *basket, int ntracks, unsigned int &error);

  int FillStackBuffer(StackLikeBuffer *buffer, int ntracks, unsigned int &error);
  
  // int AddEvent(GeantTaskData *td);
  
  /** @brief Add one event to the server */
  bool AddEvent(GeantEvent *event);

  GeantEvent *GenerateNewEvent(GeantTaskData *td, unsigned int &error);
  
  GeantEvent *ActivateEvent(GeantEvent *expected, unsigned int &error);
  
  void CompletedEvent(GeantEvent *event, GeantTaskData *td);
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_EVENT_SERVER_H
