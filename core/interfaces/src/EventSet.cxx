#include "Geant/EventSet.h"

#include <sstream>
#include "Geant/Error.h"
#include "Geant/Event.h"
#include "Geant/EventServer.h"

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

//______________________________________________________________________________
EventSet::EventSet(size_t nevents) : fNdone(0)
{
  fNevents = nevents;
  fMarkers = new EventMarker *[nevents];
  for (size_t i = 0; i < nevents; ++i) {
    fMarkers[i] = new EventMarker;
    fMarkers[i]->fDone.store(false);
  }
}

//______________________________________________________________________________
EventSet::EventSet(std::vector<Event *> const &events) : fNdone(0)
{
  fNevents = (int)events.size();
  fMarkers = new EventMarker *[fNevents];
  for (size_t i = 0; i < fNevents; ++i) {
    fMarkers[i]         = new EventMarker;
    fMarkers[i]->fEvent = events[i];
    fMarkers[i]->fDone.store(false);
  }
  fNadded = fNevents;
}

//______________________________________________________________________________
void EventSet::AddSetToServer(EventServer *evserv) const
{
  for (size_t i = 0; i < fNevents; ++i) {
    evserv->AddEvent(fMarkers[i]->fEvent);
    // Now the event number is known
    fMarkers[i]->fEventNumber = fMarkers[i]->fEvent->GetEvent();
  }
}

//______________________________________________________________________________
bool EventSet::AddEvent(Event *event)
{
  if (fNadded == fNevents) {
    Error("EventSet::AddEvent", "The event set already complete");
    return false;
  }

  fMarkers[fNadded++]->fEvent = event;
  return true;
}

//______________________________________________________________________________
void EventSet::Print() const
{
  // Print the event set content
  std::stringstream os;
  os << "events: ";
  for (size_t i = 0; i < fNevents; ++i)
    os << " " << fMarkers[i]->fEventNumber;
  std::cout << os.str() << std::endl;
}

//______________________________________________________________________________
void EventSet::Info(std::string &info) const
{
  std::stringstream os;
  for (size_t i = 0; i < fNevents; ++i) {
    os << fMarkers[i]->fEventNumber;
    if (fMarkers[i]->fDone.load())
      os << "(done)  ";
    else
      os << "  ";
  }
  info = os.str();
}

//______________________________________________________________________________
bool EventSet::MarkDone(int event_number)
{
  size_t slot;
  if (!Contains(event_number, slot)) return false;
  // Check if event is already marked as done
  if (!fMarkers[slot]->fDone.load()) {
    fMarkers[slot]->fDone.store(true);
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

} // namespace GEANT_IMPL_NAMESPACE
} // namespace geant
