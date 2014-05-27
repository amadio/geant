#include "GeantScheduler.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

#ifdef USE_VECGEOM_NAVIGATOR
 #include "navigation/simple_navigator.h"
 #include "base/vector3d.h"
 #include "management/geo_manager.h"
#endif
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"

ClassImp(GeantScheduler)

//______________________________________________________________________________
GeantScheduler::GeantScheduler()
                   :TObject(),
                    fNvolumes(0),
                    fNpriority(0),
                    fBasketMgr(0),
                    fGarbageCollector(0),
                    fNtracks(0)
{
// dummy
   fPriorityRange[0] = fPriorityRange[1] = -1;
#ifdef __STAT_DEBUG
   fPStat.InitArrays(GeantPropagator::Instance()->fNevents);
   fQStat.InitArrays(GeantPropagator::Instance()->fNevents);
   fTStat.InitArrays(GeantPropagator::Instance()->fNevents);
#endif   
}

//______________________________________________________________________________
GeantScheduler::~GeantScheduler()
{
// dtor.
   if (fBasketMgr) {
      for (Int_t ib=0; ib<fNvolumes; ib++) delete fBasketMgr[ib];
   }   
   delete [] fBasketMgr;
   delete [] fNtracks;
}

//______________________________________________________________________________
void GeantScheduler::AdjustBasketSize()
{
// Adjust the basket size to converge to Ntracks/2*Nthreads
   const Int_t min_size = 4;
   const Int_t max_size = 256;
   Int_t nthreads = gPropagator->fNthreads;
   Int_t nproposed;
   for (Int_t ib=0; ib<fNvolumes; ib++) {
      if (!fNtracks[ib]) continue;
      nproposed = fNtracks[ib]/(2*nthreads);
      nproposed -= nproposed%4;
      if (!nproposed) continue;
      if (nproposed<min_size)  nproposed = min_size;
      else if (nproposed>max_size) nproposed = max_size;
//      Printf("basket %s (%d tracks): crt=%d  proposed=%d", 
//       fBasketMgr[ib]->GetName(), fNtracks[ib], fBasketMgr[ib]->GetThreshold(), nproposed);
      fBasketMgr[ib]->SetThreshold(nproposed);
   }
}

//______________________________________________________________________________
void GeantScheduler::CreateBaskets()
{
// Create the array of baskets
   if (fBasketMgr) return;
   fNvolumes = gGeoManager->GetListOfVolumes()->GetEntries();
   fBasketMgr = new GeantBasketMgr*[fNvolumes];
   fNtracks = new Int_t[fNvolumes];
   memset(fNtracks,0,fNvolumes*sizeof(Int_t));
   dcqueue<GeantBasket> *feeder = WorkloadManager::Instance()->FeederQueue();
   TIter next(gGeoManager->GetListOfVolumes());
   TGeoVolume *vol;
   GeantBasketMgr *basket_mgr;
   Int_t icrt = 0;
   Int_t nperbasket = gPropagator->fNperBasket;
   while ((vol=(TGeoVolume*)next())) {
      basket_mgr = new GeantBasketMgr(this, vol, icrt);
      basket_mgr->SetThreshold(nperbasket);
      vol->SetFWExtension(basket_mgr);
      basket_mgr->SetFeederQueue(feeder);
//      Printf("basket %s: %p feeder=%p", basket_mgr->GetName(), basket_mgr, basket_mgr->GetFeederQueue());
      fBasketMgr[icrt++] = basket_mgr;
   }
}

//______________________________________________________________________________
Int_t GeantScheduler::AddTrack(GeantTrack &track)
{
// Main method to inject generated tracks. Track status is kNew here.
   TGeoVolume *vol = track.fPath->GetCurrentNode()->GetVolume();
   GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr*>(vol->GetFWExtension());
   fNtracks[basket_mgr->GetNumber()]++;
   return basket_mgr->AddTrack(track, kFALSE);
}

//______________________________________________________________________________
Int_t GeantScheduler::AddTracks(GeantBasket *output, Int_t &ntot, Int_t &nnew, Int_t &nkilled)
{
// Main re-dispatch method. Add all tracks and inject baskets if above threshold. 
// Returns the number of injected baskets.
   Int_t ninjected = 0;
   Bool_t priority = kFALSE;
   GeantTrack_v &tracks = output->GetOutputTracks();
   Int_t ntracks = tracks.GetNtracks();
   ntot += ntracks;
   GeantBasketMgr *basket_mgr = 0;
   Int_t output_id = output->GetBasketMgr()->GetNumber();
   TGeoVolume *vol = 0;
   for (Int_t itr=0; itr<ntracks; itr++) {
      // We have to collect the killed tracks
      if (tracks.fStatusV[itr]==kKilled ||
          tracks.fStatusV[itr]==kExitingSetup ||
          tracks.fPathV[itr]->IsOutside()) { 
         nkilled++; 
         gPropagator->StopTrack(tracks, itr); 
         tracks.DeleteTrack(itr);
         fNtracks[output_id]--;
         continue;
      }
      if (tracks.fStatusV[itr]!=kNew) fNtracks[output_id]--;
      else nnew++;
      tracks.fStatusV[itr] = kAlive;
      priority = kFALSE;
      if (fPriorityRange[0]>=0 &&
          tracks.fEventV[itr]>=fPriorityRange[0] &&
          tracks.fEventV[itr]<=fPriorityRange[1]) priority = kTRUE;
      vol = tracks.fPathV[itr]->GetCurrentNode()->GetVolume();    
      basket_mgr = static_cast<GeantBasketMgr*>(vol->GetFWExtension());
      fNtracks[output_id]++;
      ninjected += basket_mgr->AddTrack(tracks, itr, priority);
   } 
   tracks.Clear();
   return ninjected;
}  
   
//______________________________________________________________________________
Int_t GeantScheduler::CollectPrioritizedTracks()
{
// Send the signal to all basket managers to prioritize all pending tracks
// if any within the priority range.
   Int_t ninjected = 0;
   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++)
      ninjected += fBasketMgr[ibasket]->CollectPrioritizedTracks(fPriorityRange[0],fPriorityRange[1]);
   return ninjected;
}

//______________________________________________________________________________
Int_t GeantScheduler::FlushPriorityBaskets()
{
// Flush all non-empty priority baskets.
   Int_t ninjected = 0;
   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      ninjected += fBasketMgr[ibasket]->FlushPriorityBasket();
   }
   return ninjected;
}   

//______________________________________________________________________________
Int_t GeantScheduler::GarbageCollect()
{
// Flush all filled baskets in the work queue.
   Int_t ninjected = 0;
   for (Int_t ibasket=0; ibasket<fNvolumes; ibasket++) {
      ninjected += fBasketMgr[ibasket]->GarbageCollect();
   }
   return ninjected;
}
