#include "MyApplication.h"
#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "globals.h"

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
void MyApplication::StepManager(Int_t tid, Int_t iproc, Int_t npart, const GeantTrack_v & tracks)
{
// Application stepping manager. The thread id has to be used to manage storage
// of hits independently per thread.
   return;     // FOR NOW
   if (iproc != 1) return;     // hits for eloss process only
   MyHit *hit;
   Int_t nhits = 0;
   for (Int_t i=0; i<npart; i++) {
      hit = fFactory->NextFree(tracks.fEvslotV[i]);
      hit->fX = tracks.fXposV[i];
      hit->fY = tracks.fYposV[i];
      hit->fZ = tracks.fZposV[i];
      hit->fDe = 0.; // should be: tracks.fElossV[i];
      hit->fVolId = 0.;
      hit->fDetId = 0.;
//      if (track->path && track->path->GetCurrentNode()) {
//         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
//         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
//      }
      nhits++;
   }
//   Printf("Thread %d produced %d hits", tid, nhits);
}

//______________________________________________________________________________
void MyApplication::Digitize(Int_t /*event*/)
{
// User method to digitize a full event, which is at this stage fully transported
}
