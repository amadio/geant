#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
Int_t GeantEvent::AddTrack()
{
// Thread safe track addition
   Int_t itrack;
   fMutex.Lock();
   fNtracks++;
   itrack = fNtracks-1;
   fMutex.UnLock();
   return itrack;
}

//______________________________________________________________________________
void GeantEvent::StopTrack()
{
// Thread safe track addition
   fMutex.Lock();
   fNdone++;
   fMutex.UnLock();
}
