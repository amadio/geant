#ifndef GEANT_EVENT_SERVER_H
#define GEANT_EVENT_SERVER_H

#include <atomic>

#include <list>

#include "Geant/Typedefs.h"
#include "Geant/TaskData.h"
#include "Geant/Event.h"
#include "GeantConfig.h"
#include "Geant/priority_queue.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PrimaryGenerator;
class TaskData;
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

class EventServer {
public:
  enum ErrorType {
    kNoerr   = 0x00, // normal operation
    kCSlots  = 0x01, // contention on slots
    kCEvents = 0x02, // contention on queued events
    kDone    = 0x04  // all events dispatched
  };

  using queue_slots  = mpmc_bounded_queue<size_t>;
  using queue_events = mpmc_bounded_queue<Event *>;

private:
  int fNevents      = 0;         /** Number of events to be filled */
  int fNprimaries   = 0;         /** Number of primaries served */
  int fNactiveMax   = 0;         /** Maximum number of active events (buffer size)*/
  int fNtracksInit  = 0;         /** Initial number of tracks in the buffer. */
  int fNbasketsInit = 0;         /** Initial number of baskets to be served */
  int fInitialBsize = 0;         /** Initial basket size */
  int fMaxInit      = 0;         /** Maximum number of baskets taken by a thread in the initial phase */
  std::atomic_int fNactive;      /** Number of deployed events */
  std::atomic_int fNserved;      /** Number of baskets served */
  std::atomic_int fLastActive;   /** Last activated event */
  std::atomic_int fCurrentEvent; /** Current event being served */
  std::atomic_int fNload;        /** Last load event in the server */
  std::atomic_int fNstored;      /** Number of stored events in the server */
  std::atomic_int fNcompleted;   /** Number of completed events */
  std::atomic_flag fGenLock;     /** Generator lock */
  std::atomic<Event *> fEvent;   /** Current event being distributed */
  RunManager *fRunMgr = nullptr; /** Run manager */
  bool fEventsServed  = false;   /** All events served */
  // bool fDone = false;                  /** All events transported */
  bool fHasTracks    = false;   /** Server has tracks to dispatch */
  bool fInitialPhase = true;    /** Server in initial dispatch phase */
  Event **fEvents    = nullptr; /** Events to be dispatched */
  int fBindex        = 0;       /** Basket manager index */
  queue_slots fFreeSlots;       /** Queue of free event slots */
  queue_events fPendingEvents;  /** Queue of pending events */
  queue_events fDoneEvents;     /** Queue of transported events */

protected:
  Track *GetNextTrack(TaskData *td, unsigned int &error);

public:
  EventServer(int nactive_max, RunManager *runmgr);
  ~EventServer();

  GEANT_FORCE_INLINE
  size_t AdjustSize(size_t size) const
  {
    size_t i = 1, new_size = 2;
    while ((size >> i++) > 0)
      new_size *= 2;
    return new_size;
  }

  // Accessors
  GEANT_FORCE_INLINE
  int GetNevents() const { return fNevents; }

  GEANT_FORCE_INLINE
  int GetNcompleted() const { return fNcompleted; }

  GEANT_FORCE_INLINE
  int GetNprimaries() const { return fNprimaries; }

  GEANT_FORCE_INLINE
  int GetNbasketsInit() const { return fNbasketsInit; }

  GEANT_FORCE_INLINE
  int GetMaxInit() const { return fMaxInit; }

  GEANT_FORCE_INLINE
  int GetNactiveMax() const { return fNactiveMax; }

  GEANT_FORCE_INLINE
  int GetNactive() const { return fNactive.load(); }

  GEANT_FORCE_INLINE
  int GetCurrent() const { return fCurrentEvent.load(); }

  GEANT_FORCE_INLINE
  int GetNstored() const { return fNstored.load(); }

  GEANT_FORCE_INLINE
  int GetBsize() const { return fInitialBsize; }

  GEANT_FORCE_INLINE
  Event *GetEvent(int slot) const { return fEvents[slot]; }

  GEANT_FORCE_INLINE
  Event *FindEvent(int event) const
  {
    // trying to avoid concurrent map, but this could be smarter
    for (int i = 0; i < fNactiveMax; ++i) {
      if (fEvents[i]->GetEvent() == event) return fEvents[i];
    }
    return nullptr;
  }

  GEANT_FORCE_INLINE
  int GetBindex() const { return fBindex; }

  GEANT_FORCE_INLINE
  bool EventsServed() const { return fEventsServed; }

  GEANT_FORCE_INLINE
  bool HasTracks() const { return fHasTracks; }

  GEANT_FORCE_INLINE
  bool IsInitialPhase() const { return fInitialPhase; }

  int FillStackBuffer(StackLikeBuffer *buffer, int ntracks, TaskData *td, unsigned int &error);

  // int AddEvent(TaskData *td);

  /** @brief Add one event to the server */
  bool AddEvent(Event *event);

  Event *GenerateNewEvent(TaskData *td, unsigned int &error);

  Event *ActivateEvent(Event *expected, unsigned int &error, TaskData *td = nullptr);

  void CompletedEvent(Event *event, TaskData *td);
};

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant

#endif // GEANT_EVENT_SERVER_H
