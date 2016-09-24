#ifndef GEANT_EVENT_SERVER_H
#define GEANT_EVENT_SERVER_H

#include <atomic>

#ifndef GEANT_NVCC
#include <vector>
#else
#include "base/Vector.h"
#endif

#include "base/BitSet.h"
#include "Geant/Typedefs.h"
#include "GeantTaskData.h"
#include "GeantConfig.h"
#include "PrimaryGenerator.h"

using namespace Geant;
using namespace veccore;
//------------------------------------------------------------------------------
//
// A concurrent event server supporting:
// - Concurrent filling of tracks from input events
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

struct GeantInputEvent {
  template <class T>
#ifndef GEANT_NVCC
  using vector_t = std::vector<T>;
#else
  using vector_t = vecgeom::Vector<T>;
#endif  

  int fEvent = 0;                 /** Event number */
  int fNtracks = 0;               /** Number of tracks in the input event */
  vector_t<GeantTrack> fTracks;   /** Vector containing all primary tracks */
  std::atomic_int fNfilled = 0;   /** Number of tracks copied in buffer */
  std::atomic_int fNdispatched = 0; /** Number of tracks dispatched */
};
  
class GeantEventServer
{
public:
  template <class T>
#ifndef GEANT_NVCC
  using vector_t = std::vector<T>;
#else
  using vector_t = vecgeom::Vector<T>;
#endif  

private:
  int fNevents;                   /** Number of events to be filled */
  vector_t<GeantInputEvent *> fEvents; /** Events to be dispatched */
  std::atomic_int fCurrentEvent = 0;   /** Current event being served */
  std::atomic_int fLastEvent = 0;      /** Last event in the server */
  GeantRunManager *fRunMgr = nullptr;  /** Run manager */
  PrimaryGenerator *fPrimaryGenerator = nullptr; /** Primary generator */

  VolumePath_t *fStartPath;       /** Starting geometry path */

private:
  bool LoadVecGeomGeometry();
  void InitNavigators();

public:
  GeantEventPool(int event_capacity) : fNevents(event_capacity);
  ~GeantEventPool();

// Accessors
  GEANT_FORCE_INLINE
  int  GetNevents() { return fNevents; }

  GEANT_FORCE_INLINE
  void SetPrimaryGenerator(PrimaryGenerator *gen) { fPrimaryGenerator = gen; }

  bool ReadNextEvent();
  
  void AddEvent();

  void AddTrack(int event, int itrack);

  GeantTrack *GetNextTrack(int event);
};

#endif // GEANT_EVENT_SERVER_H
