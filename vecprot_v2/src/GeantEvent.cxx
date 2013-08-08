#include "GeantEvent.h"

ClassImp(GeantEvent)

//______________________________________________________________________________
void GeantEvent::AddTrack()
{
// Thread safe track addition
   the_mutex.Lock();
   ntracks++;
   the_mutex.UnLock();
}

//______________________________________________________________________________
void GeantEvent::StopTrack()
{
// Thread safe track addition
   the_mutex.Lock();
   ndone++;
   the_mutex.UnLock();
}
