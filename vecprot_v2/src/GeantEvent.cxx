#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
Int_t GeantEvent::AddTrack() {
  // Thread safe track addition
#if __cplusplus >= 201103L
  return fNtracks++;
#endif
}

//______________________________________________________________________________
void GeantEvent::Print(Option_t *) const {
  // Print events content
  Printf("Event %d: %d/%d tracks transported", GetEvent(), GetNdone(), GetNtracks());
}
