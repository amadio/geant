#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
Int_t GeantEvent::AddTrack()
{
// Thread safe track addition
   Int_t itrack;
   fMutex.Lock();
   fNtracks++;
//   itrack = fNtracks-1;
   itrack = fNtracks.load()-1;
   fMutex.UnLock();
   return itrack;
}

#if __cplusplus < 201103L
//______________________________________________________________________________
void GeantEvent::StopTrack()
{
// Thread safe track addition
   fMutex.Lock();
   fNdone++;
   fMutex.UnLock();
}
#endif

//______________________________________________________________________________
void GeantEvent::Print(Option_t *) const
{
// Print events content
  // Printf("Event %d: %d/%d tracks transported", fEvent, fNdone, fNtracks);
}
   
