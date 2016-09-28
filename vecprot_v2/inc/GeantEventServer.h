#ifndef GEANT_EVENT_SERVER_H
#define GEANT_EVENT_SERVER_H

#include <atomic>

#include <vector>

#include "base/BitSet.h"
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantConfig.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class PrimaryGenerator;
class GeantTaskData;
class GeantEvent;

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
  std::vector<GeantEvent *> fEvents; /** Events to be dispatched */
  std::atomic_int fCurrentEvent;       /** Current event being served */
  std::atomic_int fNload ;             /** Last load event in the server */
  std::atomic_int fNstored;            /** Number of stored events in the server */
  GeantRunManager *fRunMgr = nullptr;  /** Run manager */
  bool fDone = false;
 
  GeantTrack *GetNextTrack();

public:
  GeantEventServer(int event_capacity, GeantRunManager *runmgr);
  ~GeantEventServer();

// Accessors
  GEANT_FORCE_INLINE
  int  GetNevents() { return fNevents; }

  int FillBasket(GeantTrack_v &tracks, int ntracks);
  
  int AddEvent(GeantTaskData *td = nullptr);
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif // GEANT_EVENT_SERVER_H
