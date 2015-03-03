#include "GeantEvent.h"
#include <iostream>

//______________________________________________________________________________
int GeantEvent::AddTrack() {
  // Thread safe track addition
  int ntracks = ++fNtracks;
  int ninflight = ntracks - fNdone.load();
  int nmax = fNmax.load();
  // Thread safe update of the max number of tracks in flight
  while ( fNmax < ninflight && !fNmax.compare_exchange_weak(nmax, ninflight) )
    ;
  return ntracks;  
}

//______________________________________________________________________________
void GeantEvent::StopTrack() {
  // Mark one track as stopped. Check if event has to be prioritized.
  fNdone++;
  if (!fPrioritize) {
    if (GetNinflight() < fPriorityThr*GetNmax()) {
      fPrioritize = true;
      std::cout << "### Event " << GetEvent() << " prioritized at " <<
        100.*fPriorityThr << " % threshold" << std::endl;
    }  
  }
}

//______________________________________________________________________________
void GeantEvent::Print(const char *) const {
  // Print events content
  std::cout << "Event " << GetEvent() << ": " << GetNtracks() << 
    " tracks transported, max in flight " <<  GetNmax() << std::endl;
}
