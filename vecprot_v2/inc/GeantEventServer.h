#ifndef GEANT_EVENT_SERVER_H
#define GEANT_EVENT_SERVER_H

#include <atomic>

#include <list>

#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantConfig.h"

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
private:
  int fNevents = 0;                    /** Number of events to be filled */
  int fNactiveMax = 0;                 /** Maximum number of active events */
  int fNbasketsInit = 0;               /** Initial number of baskets to be served */
  std::atomic_int fNactive;            /** Number of deployed events */
  std::atomic_int fNserved;            /** Number of baskets served */
  std::atomic_int fLastActive;         /** Last activated event */
  std::atomic_int fCurrentEvent;       /** Current event being served */
  std::atomic_int fNload;              /** Last load event in the server */
  std::atomic_int fNstored;            /** Number of stored events in the server */
  std::atomic_int fNcompleted;         /** Number of completed events */
  GeantRunManager *fRunMgr = nullptr;  /** Run manager */
  bool fEventsServed = false;          /** All events served */
  bool fDone = false;                  /** All events transported */
  bool fHasTracks = false;             /** Server has tracks to dispatch */
  bool fInitialPhase = true;           /** Server in initial dispatch phase */
  GeantEvent** fEvents = nullptr;      /** Events to be dispatched */
  int  fBindex = 0;                    /** Basket manager index */

protected:
  GeantTrack *GetNextTrack();

public:
  GeantEventServer(int event_capacity, GeantRunManager *runmgr);
  ~GeantEventServer();

// Accessors
  GEANT_FORCE_INLINE
  int  GetNevents() const { return fNevents; }

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
  GeantEvent *GetEvent(int i) { return fEvents[i]; }

  GEANT_FORCE_INLINE
  int GetBindex() { return fBindex; }

  GEANT_FORCE_INLINE
  bool HasTracks() const { return fHasTracks; }

  GEANT_FORCE_INLINE
  bool IsInitialPhase() const { return fInitialPhase; }

  int FillBasket(GeantTrack_v &tracks, int ntracks);

  int FillBasket(Basket *basket, int ntracks);

  int FillStackBuffer(StackLikeBuffer *buffer, int ntracks);
  
  int AddEvent(GeantTaskData *td = nullptr);
  
  int ActivateEvents();
  
  void CompletedEvent(int evt);
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_EVENT_SERVER_H
