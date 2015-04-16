#include "TThread.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"

#include "TThread.h"
#include "TArrayI.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"

ClassImp(GeantBasket)

//______________________________________________________________________________
GeantBasket::GeantBasket()
    : TObject(), fManager(0), fTracksIn(), fTracksOut(), fAddingOp(0), fThreshold(0) {
  // Dummy constructor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
    : TObject(), fManager(mgr), fTracksIn(size, GeantPropagator::Instance()->fMaxDepth),
      fTracksOut(size, GeantPropagator::Instance()->fMaxDepth), fAddingOp(0), fThreshold(size) {
  // Default constructor.
  if (!mgr->GetVolume() || mgr->IsCollector()) 
    SetMixed(true);
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, Int_t depth)
    : TObject(), fManager(0), fTracksIn(size, depth), fTracksOut(size, depth), 
      fAddingOp(0), fThreshold(size) {
  // Default constructor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket() {
  // Destructor.
}

//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack &track) {
  // Add a new track to this basket. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(track);
#ifdef GeantV_SCH_DEBUG  
  assert(fAddingOp.load()>0);
#endif
  fTracksIn.AddTrackSync(track);
}

//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr) {
  // Add track from a track_v array. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(tracks, itr, kTRUE);
#ifdef GeantV_SCH_DEBUG  
  assert(fAddingOp.load()>0);
#endif
  fTracksIn.AddTrackSync(tracks, itr);
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend) {
  // Add multiple tracks from a track_v array.
#ifdef GeantV_SCH_DEBUG  
  assert(fAddingOp.load()>0);
#endif
  fTracksIn.AddTracks(tracks, istart, iend, kTRUE);
}

//______________________________________________________________________________
void GeantBasket::Clear(Option_t *option) {
  // Clear basket content.
  SetMixed(fManager->IsCollector());
  fTracksIn.Clear(option);
  fTracksOut.Clear(option);
}

//______________________________________________________________________________
Bool_t GeantBasket::Contains(Int_t evstart, Int_t nevents) const {
  // Checks if any of the tracks in the input array belongs to the given event 
  // range.
  return fTracksIn.Contains(evstart, nevents);
}

//______________________________________________________________________________
TGeoVolume *GeantBasket::GetVolume() const {
  // Returns volume for this basket
  return fManager->GetVolume();
}

//______________________________________________________________________________
void GeantBasket::Print(Option_t *) const {
  // Print basket content.
  Printf("*** basket %s: ninput=%3d   noutput=%3d", GetName(), GetNinput(), GetNoutput());
}

//______________________________________________________________________________
void GeantBasket::PrintTrack(Int_t /*itr*/, Bool_t /*input*/) const {
  // Print a given track.
}

//______________________________________________________________________________
void GeantBasket::Recycle() {
  // Recycle the basket to the volume scheduler.
  fManager->RecycleBasket(this);
}

//______________________________________________________________________________
void GeantBasket::SetThreshold(Int_t threshold) {
  // Set transport threshold for the basket
  if (threshold > fThreshold) {
    if (fTracksIn.Capacity() < threshold)
      fTracksIn.Resize(threshold);
    if (fTracksOut.Capacity() < threshold)
      fTracksOut.Resize(threshold);
  }
  fThreshold = threshold;
}

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

ClassImp(GeantBasketMgr)

//______________________________________________________________________________
GeantBasketMgr::GeantBasketMgr(GeantScheduler *sch, TGeoVolume *vol, Int_t number, Bool_t collector)
    : TGeoExtension(), fScheduler(sch), fVolume(vol), fNumber(number), fBcap(0), fQcap(32),
      fActive(kFALSE), fCollector(collector), fThreshold(0), fNbaskets(0), fNused(0), 
      fCBasket(0), fLock(), fQLock(), fBaskets(0), fFeeder(0), fMutex() {
  // Constructor
  fBcap = GeantPropagator::Instance()->fMaxPerBasket + 1;
  fBaskets = new mpmc_bounded_queue<GeantBasket *>(fQcap);
  // The line below to be removed when the automatic activation schema in place
  if (collector) 
    Activate();
}

//______________________________________________________________________________
GeantBasketMgr::~GeantBasketMgr() {
  // Clean up
  delete GetCBasket();
  GeantBasket *basket;
  while (fBaskets->dequeue(basket))
    delete basket;
}

//______________________________________________________________________________
void GeantBasketMgr::Activate()
{
  // Activate the manager for generating baskets.
  if (fActive) return;
  GeantBasket *basket;
  basket = new GeantBasket(fBcap, this);
  SetCBasket(basket);
  if (fCollector) {
    basket->SetMixed(true);
    Printf("Created collector basket manager");
  }  
  fActive = true;
}   

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::StealAndReplace(atomic_basket &current) {
  // Steal the current pointer content and replace with a new basket from the
  // pool. If the operation succeeds, returns the released basket which is now
  // thread local
  GeantBasket *newb, *oldb;
  // Backup cbasket. If a steal happened it does not matter
  newb = GetNextBasket(); // backup new basket to avoid leak
  oldb = current.load();  // both local vars, so identical content
  // Check if fPBasket pointer was stolen by other thread
  while (fLock.test_and_set(std::memory_order_acquire))
    ;
  Bool_t stolen = current.compare_exchange_strong(oldb, newb);
  fLock.clear(std::memory_order_release);
  if (stolen) {
    // No steal: current points now to newbasket which takes all new
    // tracks. We can push the backed-up old basket to the queue, but
    // only AFTER any ongoing track addition finishes.
    while (oldb->IsAddingOp()) {
    }; // new AddTrack go all to newbasket
    // Someone may still have a copy of the old basket  and haven't
    // started yet adding tracks to it
    // We are finally the owner of oldb
    //      assert (!oldb->IsAddingOp());
    return oldb;
  } else {
    // Current basket stolen by other thread: we need to recicle the new one
    // and ignore the backed-up old basket which was injected by the other thread
    //      assert(!newb->IsAddingOp());
    newb->Recycle();
  }
  return 0;
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::StealAndPin(atomic_basket &current) {
  // The method pins non-atomically the basket to the caller while adding the
  // fAddingOp flag. It makes sure that no StealAndReplace operation happened
  // on the pinned basket.
  GeantBasket *oldb = 0;
  while (1) {
    // the 2 lines below should be atomically coupled
//    while (fLock.test_and_set(std::memory_order_acquire))
      ;
    oldb = current.load();
    oldb->LockAddingOp();
//    fLock.clear(std::memory_order_release);
    // If no steal happened in between return basket pointer
    if (oldb == current.load())
      return oldb;
    oldb->UnLockAddingOp();
  }
}

//______________________________________________________________________________
Bool_t GeantBasketMgr::StealMatching(atomic_basket &global, GeantBasket *content) {
  // Steal the global basket if it has the right matching content
  // Prepare replacement
  GeantBasket *newb = GetNextBasket();
  if (global.compare_exchange_strong(content, newb)) {
    // Basket stolen
    while (content->IsAddingOp()) {
    };
    return kTRUE;
  } else {
    newb->Recycle();
  }
  return kFALSE;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority) {
  // Copy directly from a track_v a track to the basket manager.
  // Has to work concurrently
  // Atomically pin the basket for the adding operation
  GeantBasket *oldb = StealAndPin(fCBasket);
  // Now basket matches fP(C)Basket content and has the adding flag set
  oldb->AddTrack(trackv, itr);
  if (oldb->GetNinput() >= oldb->GetThreshold()) {
    oldb->UnLockAddingOp();
    if (StealMatching(fCBasket, oldb)) {
      // we fully own now oldb
      //         assert(!oldb->IsAddingOp());
      fFeeder->push(oldb, priority);
      return 1;
    } else
      return 0;
  }
  oldb->UnLockAddingOp();
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrackSingleThread(GeantTrack_v &trackv, Int_t itr, Bool_t priority) {
  // Copy directly from a track_v a track to the basket manager. It is 
  // assumed that this manager is only handled by a single thread.
  GeantBasket *cbasket = GetCBasket();
  cbasket->GetInputTracks().AddTrack(trackv, itr);
  if (cbasket->GetNinput() >= cbasket->GetThreshold()) {
    fFeeder->push(cbasket, priority);    
    SetCBasket(GetNextBasket());
    return 1;
  }
  return 0;
}
  
//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack &track, Bool_t priority) {
  // Add a track to the volume basket manager. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. The feeder must be defined beforehand. Returns the number of dispatched
  // baskets
  // Has to work concurrently
  // Atomically pin the basket for the adding operation
  GeantBasket *oldb = StealAndPin(fCBasket);
  // Now basket matches fP(C)Basket content and has the adding flag set
  oldb->AddTrack(track);
  if (oldb->GetNinput() >= oldb->GetThreshold()) {
    oldb->UnLockAddingOp();
    if (StealMatching(fCBasket, oldb)) {
      // we fully own now oldb
      //         assert(!oldb->IsAddingOp());
      fFeeder->push(oldb, priority);
      return 1;
    } else
      return 0;
  }
  oldb->UnLockAddingOp();
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracksNew(GeantBasketMgr *gc) {
  // Garbage collect tracks from the given event range. The basket gc should
  // be thread local
  if (!GetCBasket()->GetNinput())
    return 0;
  GeantBasket *cbasket = 0;
  GeantPropagator *propagator = GeantPropagator::Instance();
  while (cbasket == 0)
    cbasket = StealAndReplace(fCBasket);
  // cbasket is now thread local
  // Loop all tracks in the current basket and mark for copy the ones belonging
  // to the desired events
  GeantTrack_v &tracks = cbasket->GetInputTracks();
  Int_t ntracks = tracks.GetNtracks();
  Int_t icoll = 0;
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (propagator->fEvents[tracks.fEvslotV[itr]]->IsPrioritized() ) {
      tracks.MarkRemoved(itr);
      icoll++;
    }
  }
  if (icoll) tracks.Compact(&gc->GetCBasket()->GetInputTracks());
  return icoll;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracks(Int_t evmin, Int_t evmax) {
  // Move current basket tracks to priority one.
  // *** NONE *** This should be done for all basket managers only once when
  // starting prioritizing an event range.
  // Lock and swap containers
  if (!GetCBasket()->GetNinput())
    return 0;
  //   Printf("=== CollectPrioritized");
  //   Print();
  GeantBasket *cbasket = 0, *basket;
  Int_t npush = 0;
  // We want to steal fCBasket
  while (cbasket == 0)
    cbasket = StealAndReplace(fCBasket);
  // cbasket is now thread local
  GeantTrack_v &tracks = cbasket->GetInputTracks();
  Int_t ntracks = tracks.GetNtracks();
  // Inject cbasket if it contains tracks with priority 
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (tracks.fEventV[itr] >= evmin && tracks.fEventV[itr] <= evmax) {
      fFeeder->push(cbasket, kTRUE);
      return 1;
    }
  }
  // if cbasket empty -> just recycle
  if (cbasket->GetNinput() == 0) {
    cbasket->Recycle();
    //         Print();
    return npush;
  }
  // cbasket has to be merged back into fCBasket. Most profitable is to
  // exchange the pointer held by fCBasket (which is most likely empty if
  // no other thread was adding on top) with cbasket.
  cbasket = fCBasket.exchange(cbasket);
  while (cbasket->IsAddingOp()) {
  }; // wait possible adding to finish
  // Now add cbasket content on top of fCBasket
  GeantTrack_v &ctracks = cbasket->GetInputTracks();
  ntracks = ctracks.GetNtracks();
  for (Int_t itr = 0; itr < ntracks; itr++) {
    basket = StealAndPin(fCBasket);
    basket->AddTrack(ctracks, itr);
    basket->UnLockAddingOp();
  }
  ctracks.Clear();
  cbasket->Recycle();
  //   Print();
  return npush;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::FlushPriorityBasket() {
  // Flush the baskets containing tracks. Returns the number of dispatched baskets.
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::GarbageCollect() {
  // Copy all priority tracks to the current basket and flush to queue
  GeantBasket *cbasket = 0;
  // We want to steal fCBasket
  if (GetCBasket()->GetNinput()) {
    while (cbasket == 0)
      cbasket = StealAndReplace(fCBasket);
    Int_t ntracks = cbasket->GetNinput();
    if (ntracks) {
      //            assert(!cbasket->IsAddingOp());
      fFeeder->push(cbasket, kFALSE);
      return 1;
    }
    cbasket->Recycle();
  }
  return 0;
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket() {
  // Returns next empy basket if any available, else create a new basket.
  //   GeantBasket *next = fBaskets->try_pop();
  GeantBasket *next = 0;
  const Int_t nthreads = GeantPropagator::Instance()->fNthreads;
  Int_t threshold = fThreshold.load();
  Int_t threshold_new = threshold * fNused.load() / nthreads;
  if ((!fCollector) && (threshold_new < fBcap) && ((threshold_new - threshold) > 4)) {
    Int_t remainder = threshold_new % 4;
    if (remainder > 0)
      threshold_new += 4 - remainder;
    fThreshold.store(threshold_new);
    //      Printf("volume %s:      ntotal=%d   nused=%d   thres=%d", GetName(), fNbaskets.load(),
    //      fNused.load(), threshold_new);
  }
  Bool_t pulled = fBaskets->dequeue(next);
  if (!pulled) {
    next = new GeantBasket(fBcap, this);
    if (fCollector) next->SetMixed(kTRUE);
    // === critical section if atomics not supported ===
    fNbaskets++;
    fNused++;
    // === end critical section ===
  } else {
    // === critical section if atomics not supported ===
    fNused++;
    // === end critical section ===
  }
  next->SetThreshold(fThreshold.load());
  return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b) {
  // Recycle a basket.
  //   assert(!b->GetNinput());
  //   assert(!b->IsAddingOp());
  const Int_t nthreads = GeantPropagator::Instance()->fNthreads;
  Int_t threshold = fThreshold.load();
  Int_t threshold_new = threshold * fNused.load() / nthreads;
  if ((!fCollector) && (threshold_new > 4) && ((threshold_new - threshold) < -4)) {
    Int_t remainder = threshold_new % 4;
    threshold_new -= remainder;
    fThreshold.store(threshold_new);
    //      Printf("volume %s:      ntotal=%d   nused=%d   thres=%d    size=%ld", GetName(),
    //      fNbaskets.load(), fNused.load(), threshold_new, Sizeof());
  }
  b->Clear();
  Int_t nbaskets = fNbaskets.load();
  //   Int_t nused = fNused.load();
  if (nbaskets > fQcap) {
    //      Printf("######## Increasing queue size for %s to %d ########", GetName(), 2*fQcap);
    IncreaseQueueSize(2 * fQcap);
  }
  if (!fBaskets->enqueue(b)) {
    // The queue is full - > delete the basket
    //      Printf("=== ALARM %s: threshold=%d", GetName(), fThreshold.load());
    //      Printf("Deleting basket %p for %s", (void*)b, GetName());
    delete b;
    fNbaskets--;
    //      Printf("Fatal error: exceeded the size of the bounded queue for basket mgr: %s",
    //      b->GetName());
    //      exit(1);
  }
  fNused--;
}
//______________________________________________________________________________
void GeantBasketMgr::CleanBaskets(Int_t ntoclean) {
  // Clean a number of recycled baskets to free some memory
  GeantBasket *toclean;
  for (auto i = 0; i < ntoclean; i++) {
    if (fBaskets->dequeue(toclean)) {
      delete toclean;
      fNbaskets--;
    }
  }
}

//______________________________________________________________________________
Bool_t GeantBasketMgr::IncreaseQueueSize(Int_t newsize) {
  // Increase dynamically the queue size
  auto oldbaskets = fBaskets;
  if (fQLock.test_and_set(std::memory_order_acquire))
    return kFALSE;
  if (newsize <= fQcap)
    return kFALSE;
  fBaskets = new mpmc_bounded_queue<GeantBasket *>(newsize);
  fQcap = newsize;
  fQLock.clear(std::memory_order_release);
  GeantBasket *next;
  while (oldbaskets->dequeue(next))
    if (!fBaskets->enqueue(next))
      delete next;
  return kTRUE;
}

//______________________________________________________________________________
void GeantBasketMgr::Print(Option_t *) const {
  // Print info about the basket content.
  Printf("Bsk_mgr %s: current: in=%d out=%d", GetName(),
         GetCBasket()->GetNinput(), GetCBasket()->GetNoutput());
}

//______________________________________________________________________________
void GeantBasketMgr::PrintSize() const {
  // Print detailed info about size.
  size_t size = Sizeof();
  size_t sizeb = 0;
  if (GetCBasket())
    sizeb = GetCBasket()->Sizeof();
  Printf("Bsk_mgr %s: %d baskets of size %ld:    %ld bytes", GetName(), GetNbaskets(), sizeb, size);
}
