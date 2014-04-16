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
               fInitialized(kFALSE),
               fIdGap(0),
               fIdAbs(0),
               fFactory(0)
{
// Ctor..
   GeantFactoryStore *store = GeantFactoryStore::Instance();
   fFactory = store->GetFactory<MyHit>(16); 
   printf("Created factory for MyHit");
   memset(fEdepGap, 0, kNlayers*sizeof(Float_t));
   memset(fLengthGap, 0, kNlayers*sizeof(Float_t));
   memset(fEdepAbs, 0, kNlayers*sizeof(Float_t));
   memset(fLengthAbs, 0, kNlayers*sizeof(Float_t));
}

//______________________________________________________________________________
Bool_t MyApplication::Initialize()
{
// Initialize application. Geometry must be loaded.
   if (fInitialized) return kTRUE;
   if (!gGeoManager) {
      Error("Initialize", "Geometry not loaded");
      return kFALSE;
   }
   TGeoVolume *lvGap = gGeoManager->GetVolume("liquidArgon");
   TGeoVolume *lvAbs = gGeoManager->GetVolume("Lead");
   
   if (!lvGap || !lvAbs) {
      Error("Initialize","Logical volumes for gap and absorber not found - do you use the right geometry");
      return kFALSE;
   }
   fIdGap = lvGap->GetNumber();
   fIdAbs = lvAbs->GetNumber();
   fInitialized = kTRUE;
   return kTRUE;
}   
   
//______________________________________________________________________________
void MyApplication::StepManager(Int_t tid, Int_t npart, const GeantTrack_v & tracks)
{
// Application stepping manager. The thread id has to be used to manage storage
// of hits independently per thread.
   if (!fInitialized) return;     // FOR NOW
   // Loop all tracks, check if they are in the right volume and collect the
   // energy deposit and step length
   TGeoNode *current;
   Int_t idvol, idnode, ilev;
   for (Int_t i=0; i<npart; i++) {
//      printf("%d=>\n", i);
//      tracks.PrintTrack(i);
      ilev = tracks.fPathV[i]->GetLevel();
      if (ilev<1) continue;
      current = tracks.fPathV[i]->GetCurrentNode();
      if (!current) continue;
      idnode = tracks.fPathV[i]->GetNode(ilev-1)->GetNumber();
      idvol = current->GetVolume()->GetNumber();
      if (idvol == fIdGap) {
//         tracks.PrintTrack(i);
         fEdepGap[idnode-3] += tracks.fEdepV[i];
         fLengthGap[idnode-3] += tracks.fStepV[i];
      } else if (idvol == fIdAbs) {
//         tracks.PrintTrack(i);
         fEdepAbs[idnode-3] += tracks.fEdepV[i];
         fLengthAbs[idnode-3] += tracks.fStepV[i];
      }
   }   
   return;
   MyHit *hit;
   Int_t nhits = 0;
   for (Int_t i=0; i<npart; i++) {
      hit = fFactory->NextFree(tracks.fEvslotV[i]);
      hit->fX = tracks.fXposV[i];
      hit->fY = tracks.fYposV[i];
      hit->fZ = tracks.fZposV[i];
      hit->fEdep = tracks.fEdepV[i];
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
void MyApplication::Digitize(Int_t event)
{
// User method to digitize a full event, which is at this stage fully transported
   printf("======= Statistics for event %d:\n", event);
   for (Int_t i=0; i<kNlayers; ++i) {
      printf("Layer %d: Egap=%f   Lgap=%f   Eabs=%f   Labs=%f\n", i+3, 
             fEdepGap[i], fLengthGap[i], fEdepAbs[i], fLengthAbs[i]);
   }
//   memset(fEdepGap, 0, kNlayers*sizeof(Float_t));
//   memset(fLengthGap, 0, kNlayers*sizeof(Float_t));
//   memset(fEdepAbs, 0, kNlayers*sizeof(Float_t));
//   memset(fLengthAbs, 0, kNlayers*sizeof(Float_t));
}
