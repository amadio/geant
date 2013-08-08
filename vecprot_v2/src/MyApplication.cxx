#include "MyApplication.h"
#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"

ClassImp(MyApplication)

//______________________________________________________________________________
MyApplication::MyApplication()
              :GeantVApplication(),
               fFactory(0)
{
// Ctor..
   GeantFactoryStore *store = GeantFactoryStore::Instance();
   fFactory = store->GetFactory<MyHit>(16); 
}

//______________________________________________________________________________
void MyApplication::StepManager(Int_t tid, Int_t iproc, Int_t npart, Int_t *particles)
{
// Application stepping manager. The thread id has to be used to manage storage
// of hits independently per thread.
   if (iproc != 1) return;     // hits for eloss process only
   GeantTrack *track;
   MyHit *hit;
   Int_t nhits = 0;
   for (Int_t i=0; i<npart; i++) {
      track = gPropagator->fTracks[particles[i]];
      hit = fFactory->NextFree(track->evslot);
      hit->fX = track->xpos;
      hit->fY = track->ypos;
      hit->fZ = track->zpos;
      hit->fDe = track->eloss;
      if (track->path && track->path->GetCurrentNode()) {
         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
      }
      nhits++;
   }
//   Printf("Thread %d produced %d hits", tid, nhits);
}

//______________________________________________________________________________
void MyApplication::Digitize(Int_t /*event*/)
{
// User method to digitize a full event, which is at this stage fully transported
}
