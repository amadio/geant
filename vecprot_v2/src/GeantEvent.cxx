#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
void GeantEvent::AddTrack()
{
// Thread safe track addition
   fMutex.Lock();
   fNtracks++;
   fMutex.UnLock();
}

//______________________________________________________________________________
void GeantEvent::StopTrack()
{
// Thread safe track addition
   fMutex.Lock();
   fNdone++;
   fMutex.UnLock();
}
