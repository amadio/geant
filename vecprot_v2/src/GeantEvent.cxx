#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
Int_t GeantEvent::AddTrack() {
  // Thread safe track addition
  Int_t ntracks = ++fNtracks;
  Int_t ninflight = ntracks - fNdone.load();
  Int_t nmax = fNmax.load();
  // Thread safe update of the max number of tracks in flight
  while ( fNmax < ninflight && !fNmax.compare_exchange_weak(nmax, ninflight) )
    ;
  return ntracks;  
}

//______________________________________________________________________________
void GeantEvent::Print(Option_t *) const {
  // Print events content
  Printf("Event %d: %d tracks transported, %d max in flight", GetEvent(), GetNtracks(), GetNmax());
}
