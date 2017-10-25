#include "EventSet.h"

#include "Geant/Error.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
EventSet::EventSet(size_t nevents) : fNdone(0) {
  fNevents = nevents;
  fMarkers = new EventMarker*[nevents];
  for (size_t i = 0; i < nevents; ++i) {
    fMarkers[i] = new EventMarker;
    fMarkers[i]->fEvent = -1;
    fMarkers[i]->fDone.clear();
  }
}

//______________________________________________________________________________
EventSet::EventSet(std::vector<int> const &events) : fNdone(0) {
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
bool EventSet::AddEvent(int event) {
  if (fNadded == fNevents) {
    Error("EventSet::AddEvent", "The event set already complete");
    return false;
  }
  
  size_t slot;
  if (Contains(event, slot)) {
    Error("EventSet::AddEvent", "The event set already contains event %d", event);
    return false;  
  }
  
  fMarkers[fNadded++]->fEvent = event;
  return true;
}

//______________________________________________________________________________
bool EventSet::MarkDone(int event) {
  size_t slot;
  if (!Contains(event, slot)) return false;
  // Check if event is already marked as done
  if (!fMarkers[slot]->fDone.test_and_set(std::memory_order_acquire)) {
    size_t ndone = fNdone.fetch_add(1) + 1;
    if (ndone == fNevents) {
      fDone = true;
      return true;
    }
  }
  return false;
}

} // GEANT_IMPL_NAMESPACE
} // Geant
