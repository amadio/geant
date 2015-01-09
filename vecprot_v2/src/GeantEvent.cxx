#include "GeantEvent.h"

ClassImp(GeantEvent)

    //______________________________________________________________________________
    Int_t GeantEvent::AddTrack() {
// Thread safe track addition
#if __cplusplus >= 201103L
  return fNtracks++;
#else
  fMutex.Lock();
  Int_t itrack = fNtracks++;
  fMutex.UnLock();
  return itrack;
#endif
}

//______________________________________________________________________________
void GeantEvent::Print(Option_t *) const {
  // Print events content
  // Printf("Event %d: %d/%d tracks transported", fEvent, fNdone, fNtracks);
}
