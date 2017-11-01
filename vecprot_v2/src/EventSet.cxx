#include "EventSet.h"

#include "Geant/Error.h"
#include "GeantEvent.h"
#include "GeantEventServer.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
EventSet::EventSet(size_t nevents) : fNdone(0) {
  fNevents = nevents;
  fMarkers = new EventMarker*[nevents];
  for (size_t i = 0; i < nevents; ++i) {
    fMarkers[i] = new EventMarker;
    fMarkers[i]->fDone.clear();
  }
}

//______________________________________________________________________________
EventSet::EventSet(std::vector<GeantEvent*> const &events) : fNdone(0) {
  fNevents = (int)events.size();
  fMarkers = new EventMarker*[fNevents];
  for (size_t i = 0; i < fNevents; ++i) {
    fMarkers[i] = new EventMarker;
    fMarkers[i]->fEvent = events[i];
    fMarkers[i]->fDone.clear();
  }
  fNadded = fNevents;
}

//______________________________________________________________________________
void EventSet::AddSetToServer(GeantEventServer *evserv) const
{
  for (size_t i = 0; i < fNevents; ++i)
    evserv->AddEvent(fMarkers[i]->fEvent);
}

//______________________________________________________________________________
bool EventSet::AddEvent(GeantEvent *event) {
  if (fNadded == fNevents) {
    Error("EventSet::AddEvent", "The event set already complete");
    return false;
  }
  size_t slot;
  if (Contains(event, slot)) {
    Error("EventSet::AddEvent", "The event set already contains same event %d", event->GetEvent());
    return false;  
  }
  
  fMarkers[fNadded++]->fEvent = event;
  return true;
}

//______________________________________________________________________________
bool EventSet::MarkDone(GeantEvent *event) {
  size_t slot;
  if (!Contains(event, slot)) return false;
  // Check if event is already marked as done
  if (!fMarkers[slot]->fDone.test_and_set(std::memory_order_acquire)) {
    size_t ndone = fNdone.fetch_add(1) + 1;
    if (ndone == fNevents) {
      fDone = true;
      // Wake-up all threads needing this event set done
      fLocker.StartAll();
      return true;
    }
  }
  return true;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
