#include "GeantScheduler.h"

#include "Geant/Error.h"
#include "TMath.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantTaskData.h"
#include "GeantEvent.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/SimpleNavigator.h"
#include "base/Vector3D.h"
#include "management/GeoManager.h"
using vecgeom::GeoManager;
#else
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TGeoManager.h"
#endif

//______________________________________________________________________________
GeantScheduler::GeantScheduler()
    : fNvolumes(0), fNpriority(0), fBasketMgr(0), fGarbageCollector(0), fNstvol(0), fIstvol(0), fNvect(0),
      fNsteps(0), fCrtMgr(0), fCollecting(false), fLearning(ATOMIC_FLAG_INIT), fGBCLock(ATOMIC_FLAG_INIT),
      fVolumes() {
  // Default constructor
  SetLearning(false);
  fNvect = new int[257];
  memset(fNvect, 0, 257 * sizeof(int));
}

//______________________________________________________________________________
GeantScheduler::~GeantScheduler() {
  // dtor.
  delete fGarbageCollector;
  if (fBasketMgr) {
    for (int ib = 0; ib < fNvolumes; ib++) {
#ifdef USE_VECGEOM_NAVIGATOR
      fBasketMgr[ib]->GetVolume()->SetBasketManagerPtr(0);
#else
      fBasketMgr[ib]->GetVolume()->SetFWExtension(0);
#endif
      delete fBasketMgr[ib];
    }
  }
  delete[] fBasketMgr;
  delete[] fNstvol;
  delete[] fIstvol;
  delete[] fNvect;
}

//______________________________________________________________________________
void GeantScheduler::ActivateBasketManagers() {
  // Activate basket managers based on the distribution of steps in corresponding
  // volumes.
  TMath::Sort(fNvolumes, fNstvol, fIstvol);
  int ntot = 0;
  int nsum = 0;
  // Percent of steps to activate basketized volumes
  double threshold = 0.9;
  int nactive = 0;
  const Volume_t *vol;
  for (auto i = 0; i < fNvolumes; i++) {
    ntot += fNstvol[i];
    vol = fVolumes[i];
#ifdef USE_VECGEOM_NAVIGATOR
    GeantBasketMgr *mgr = (GeantBasketMgr *)vol->GetBasketManagerPtr();
#else
    GeantBasketMgr *mgr = (GeantBasketMgr *)vol->GetFWExtension();
#endif
    if (mgr->IsActive()) {
      nsum += fNstvol[i];
      nactive++;
    }
  }
  int nthreshold = ntot * threshold;
  for (auto i = 0; i < fNvolumes; ++i) {
    vol = fVolumes[fIstvol[i]];
#ifdef USE_VECGEOM_NAVIGATOR
    GeantBasketMgr *mgr = (GeantBasketMgr *)vol->GetBasketManagerPtr();
#else
    GeantBasketMgr *mgr = (GeantBasketMgr *)vol->GetFWExtension();
#endif
    if (!mgr->IsActive()) {
      mgr->Activate();
      nsum += fNstvol[fIstvol[i]];
      nactive++;
    }
    if (nsum > nthreshold)
      break;
  }
  threshold = double(nsum) / ntot;
  Geant::Info("ActivateBasketManagers", "Activated %d volumes accounting for %4.1f%% of track steps", nactive, 100 * threshold);
  int nprint = 10;
  if (nprint > fNvolumes)
    nprint = fNvolumes;
  for (auto i = 0; i < nprint; ++i) {
    vol = fVolumes[fIstvol[i]];
    Geant::Print(" *", " %s: %d steps", vol->GetName(), fNstvol[fIstvol[i]]);
  }
}

//______________________________________________________________________________
void GeantScheduler::AdjustBasketSize() {
  // Adjust the basket size to converge to Ntracks/2*Nthreads
  return; // !!!!!!!!!!!!!!! NOT needed anymore !!!!!!!!!!!!!!!
          /*
            const int min_size = 4;
            const int max_size = gPropagator->fNperBasket;
            int nthreads = gPropagator->fNthreads;
            int nproposed;
            for (int ib = 0; ib < fNvolumes; ib++) {
              if (!GetNtracks(ib))
                continue;
              nproposed = GetNtracks(ib) / (2 * nthreads);
              nproposed -= nproposed % 4;
              if (!nproposed)
                continue;
              if (nproposed < min_size)
                nproposed = min_size;
              else if (nproposed > max_size)
                nproposed = max_size;
              //      Printf("basket %s (%d tracks): crt=%d  proposed=%d",
              //       fBasketMgr[ib]->GetName(), fNtracks[ib], fBasketMgr[ib]->GetThreshold(), nproposed);
              fBasketMgr[ib]->SetThreshold(nproposed);
            }
          */
}

//______________________________________________________________________________
void GeantScheduler::CreateBaskets() {
  // Create the array of baskets
  if (fBasketMgr)
    return;
#ifdef USE_VECGEOM_NAVIGATOR
  GeoManager::Instance().GetAllLogicalVolumes(fVolumes);
  fNvolumes = fVolumes.size();
#else
  TObjArray *lvolumes = gGeoManager->GetListOfVolumes();
  fNvolumes = lvolumes->GetEntries();
  for (auto ivol=0; ivol<fNvolumes; ivol++) fVolumes.push_back((TGeoVolume*)lvolumes->At(ivol));
#endif
  fBasketMgr = new GeantBasketMgr *[fNvolumes];
  fNstvol = new int[fNvolumes];
  fIstvol = new int[fNvolumes];
  memset(fNstvol, 0, fNvolumes * sizeof(int));
  memset(fIstvol, 0, fNvolumes * sizeof(int));
  Geant::priority_queue<GeantBasket *> *feeder = WorkloadManager::Instance()->FeederQueue();
  Volume_t *vol;
  GeantBasketMgr *basket_mgr;
  int icrt = 0;
  int nperbasket = gPropagator->fNperBasket;
  for (auto ivol=0; ivol<fNvolumes; ++ivol) {
    vol = (Volume_t *)fVolumes[ivol];
    basket_mgr = new GeantBasketMgr(this, vol, icrt);
    basket_mgr->SetThreshold(nperbasket);
#ifdef USE_VECGEOM_NAVIGATOR
    vol->SetBasketManagerPtr(basket_mgr);
#else
    vol->SetFWExtension(basket_mgr);
#endif
    basket_mgr->SetFeederQueue(feeder);
    //      Printf("basket %s: %p feeder=%p", basket_mgr->GetName(), basket_mgr,
    //      basket_mgr->GetFeederQueue());
    fBasketMgr[icrt++] = basket_mgr;
  }
  //   const size_t MB = 1048576;
  // size_t size = Sizeof() / MB;
  // Printf("Size of scheduler including initial baskets: %ld MB", size);
  //   PrintSize();
}

//______________________________________________________________________________
int GeantScheduler::AddTrack(GeantTrack &track, GeantTaskData *td) {
// Main method to inject generated tracks. Track status is kNew here.
#ifdef USE_VECGEOM_NAVIGATOR
  Volume_t *vol = const_cast<Volume_t *>(track.fPath->Top()->GetLogicalVolume());
  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->GetBasketManagerPtr());
#else
  Volume_t *vol = track.fPath->GetCurrentNode()->GetVolume();
  GeantBasketMgr *basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
#endif
  int ivol = basket_mgr->GetNumber();
  fNstvol[ivol]++;
  fNsteps++;
  // If no learning phase requested, activate all basket managers
  GeantPropagator *propagator = GeantPropagator::Instance();
  if ((propagator->fLearnSteps == 0) && !fLearning.test_and_set(std::memory_order_acquire)) {
    for (ivol = 0; ivol < fNvolumes; ++ivol)
      fBasketMgr[ivol]->Activate();
  } else {
    // Activate the basket manager
    basket_mgr->Activate();
  }
  return basket_mgr->AddTrack(track, false, td);
}

//______________________________________________________________________________
int GeantScheduler::ReusableTracks(GeantTrack_v &tracks) const
{
// Check if the basket can be reused efficiently in the next transport iteration.
  int ntracks = tracks.GetNtracks();
  int nreusable = 0;
  for (int itr = 0; itr < ntracks; ++itr) {
    if (tracks.fStatusV[itr] == kKilled || tracks.fStatusV[itr] == kExitingSetup || tracks.fPathV[itr]->IsOutside())
      continue;
    if (tracks.fStatusV[itr] == kNew || tracks.fStatusV[itr] == kPhysics)
      nreusable++;
  }
  return nreusable;
}

//______________________________________________________________________________
int GeantScheduler::CopyReusableTracks(GeantTrack_v &tracks, GeantTrack_v &input, int nmax) const
{
// Copy reusable tracks from the output tracks to the input tracks, not
// exceeding nmax
  int ntracks = tracks.GetNtracks();
  int nreused = 0;
  for (int itr=0; itr<ntracks; ++itr) {
    if (tracks.fStatusV[itr] == kKilled || tracks.fStatusV[itr] == kExitingSetup || tracks.fPathV[itr]->IsOutside())
      continue;
    if (tracks.fStatusV[itr] == kNew || tracks.fStatusV[itr] == kPhysics) {
      tracks.MarkRemoved(itr);
      nreused++;
      if (nreused == nmax) break;
    }
  }
  tracks.Compact(&input);
  return nreused;
}

//______________________________________________________________________________
int GeantScheduler::AddTracks(GeantTrack_v &tracks, int &ntot, int &nnew, int &nkilled, GeantTaskData *td) {
  // Main re-basketizing method. Add all tracks and inject baskets if above threshold.
  // Returns the number of injected baskets.
  int ninjected = 0;
  bool priority = kFALSE;
  GeantPropagator *propagator = GeantPropagator::Instance();
  int ntracks = tracks.GetNtracks();
  ntot += ntracks;
  GeantBasketMgr *basket_mgr = 0;
  Volume_t *vol = 0;
  for (int itr = 0; itr < ntracks; ++itr) {
    // We have to collect the killed tracks
    if (tracks.fStatusV[itr] == kKilled || tracks.fStatusV[itr] == kExitingSetup || tracks.fPathV[itr]->IsOutside()) {
      nkilled++;
      gPropagator->StopTrack(tracks, itr);
      tracks.DeleteTrack(itr);
      continue;
    }
    if (tracks.fStatusV[itr] == kNew)
      nnew++;
    tracks.fStatusV[itr] = kAlive;

    priority = kFALSE;

#ifdef USE_VECGEOM_NAVIGATOR
    vol = const_cast<Volume_t *>(tracks.fPathV[itr]->Top()->GetLogicalVolume());
    basket_mgr = static_cast<GeantBasketMgr *>(vol->GetBasketManagerPtr());
#else
    vol = tracks.fPathV[itr]->GetCurrentNode()->GetVolume();
    basket_mgr = static_cast<GeantBasketMgr *>(vol->GetFWExtension());
#endif
    int ivol = basket_mgr->GetNumber();
//    tracks.fVindexV[itr] = ivol;
    fNstvol[ivol]++;
    long nsteps = ++fNsteps;
    // Detect if the event the track is coming from is prioritized
    if (propagator->fEvents[tracks.fEvslotV[itr]]->IsPrioritized()) {
      ninjected += td->fBmgr->AddTrackSingleThread(tracks, itr, true, td);
      continue;
    }
    if (propagator->fLearnSteps && (nsteps % propagator->fLearnSteps) == 0 &&
        !fLearning.test_and_set(std::memory_order_acquire)) {
      Geant::Info("AddTracks", "=== Learning phase of %d steps completed ===", propagator->fLearnSteps);
      // Here comes the algorithm activating basket managers...
      ActivateBasketManagers();
      propagator->fLearnSteps *= 4;
      fLearning.clear();
    }
    if (basket_mgr->IsActive())
      ninjected += basket_mgr->AddTrack(tracks, itr, priority, td);
    else
      ninjected += td->fBmgr->AddTrackSingleThread(tracks, itr, false, td);
  }
  tracks.Clear();
  return ninjected;
}

//______________________________________________________________________________
int GeantScheduler::GarbageCollect(GeantTaskData *td, bool force) {
  // Flush all filled baskets in the work queue.
  //   PrintSize();
  // Only one thread at a time
  if (force) {
    while (fGBCLock.test_and_set(std::memory_order_acquire))
      ;
  } else {
    if (fGBCLock.test_and_set(std::memory_order_acquire))
      return 0;
  }
  int ninjected = 0;
  for (int ibasket = 0; ibasket < fNvolumes; ibasket++) {
    if (fBasketMgr[ibasket]->IsActive())
      ninjected += fBasketMgr[ibasket]->GarbageCollect(td);
  }
  fGBCLock.clear(std::memory_order_release);
//  Printf("=== Garbage collect: %d baskets", ninjected);
  return ninjected;
}

//______________________________________________________________________________
size_t GeantScheduler::Sizeof() const {
  // Returns size of the scheduler, including allocated baskets
  size_t size = sizeof(GeantScheduler);
  for (auto i = 0; i < fNvolumes; ++i)
    size += fBasketMgr[i]->Sizeof();
  if (fGarbageCollector)
    size += fGarbageCollector->Sizeof();
  return size;
}

//______________________________________________________________________________
void GeantScheduler::PrintSize() const {
  // Prints detailed breakdown of size allocated
  size_t size = Sizeof();
  Geant::Print("","Size of scheduler: %ld bytes", size);
  for (auto i = 0; i < fNvolumes; ++i)
    fBasketMgr[i]->PrintSize();
  if (fGarbageCollector)
    fGarbageCollector->PrintSize();
}
