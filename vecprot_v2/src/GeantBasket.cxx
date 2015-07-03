#include "TThread.h"
#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "PhysicsProcess.h"
#include "WorkloadManager.h"

#include "TThread.h"
#include "TArrayI.h"
#ifdef USE_VECGEOM_NAVIGATOR
#else
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"
#endif

ClassImp(GeantBasket)

    //______________________________________________________________________________
    GeantBasket::GeantBasket()
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
      fThreshold(0), fTracksIn(), fTracksOut() {
  // Dummy constructor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
    : TObject(), fManager(mgr), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
      fThreshold(size), fTracksIn(size, GeantPropagator::Instance()->fMaxDepth),
      fTracksOut(size, GeantPropagator::Instance()->fMaxDepth) {
  // Default constructor.
  if (!mgr->GetVolume() || mgr->IsCollector())
    SetMixed(true);
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, Int_t depth)
<<<<<<< HEAD
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
      fThreshold(size), fTracksIn(size, depth), fTracksOut(size, depth) {
=======
    : TObject(), fManager(0), fTracksIn(size, depth), fTracksOut(size, depth), fAddingOp(0), fThreshold(size) {
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  // Default constructor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket() {
  // Destructor.
}

//______________________________________________________________________________
<<<<<<< HEAD
Int_t GeantBasket::AddTrack(GeantTrack &track) {
  // Add a new track to this basket. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(track);
  return (fTracksIn.AddTrackSync(track));
}

//______________________________________________________________________________
Int_t GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr) {
  // Add track from a track_v array. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(tracks, itr, kTRUE);
  return (fTracksIn.AddTrackSync(tracks, itr));
=======
void GeantBasket::AddTrack(GeantTrack &track) {
// Add a new track to this basket. Has to work concurrently.

// Activating the line below adds non-concurrently an input track
//   fTracksIn.AddTrack(track);
#ifdef GeantV_SCH_DEBUG
  assert(fAddingOp.load() > 0);
#endif
  fTracksIn.AddTrackSync(track);
}

//______________________________________________________________________________
void GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr) {
// Add track from a track_v array. Has to work concurrently.

// Activating the line below adds non-concurrently an input track
//   fTracksIn.AddTrack(tracks, itr, kTRUE);
#ifdef GeantV_SCH_DEBUG
  assert(fAddingOp.load() > 0);
#endif
  fTracksIn.AddTrackSync(tracks, itr);
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend) {
<<<<<<< HEAD
  // Add multiple tracks from a track_v array.
=======
// Add multiple tracks from a track_v array.
#ifdef GeantV_SCH_DEBUG
  assert(fAddingOp.load() > 0);
#endif
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  fTracksIn.AddTracks(tracks, istart, iend, kTRUE);
}

//______________________________________________________________________________
void GeantBasket::Clear(Option_t *option) {
  // Clear basket content.
  SetMixed(fManager->IsCollector());
  fTracksIn.Clear(option);
  fTracksOut.Clear(option);
  fNbooked.store(0, std::memory_order_relaxed);
  fNcopying.store(0, std::memory_order_relaxed);
  fNcopied.store(0, std::memory_order_relaxed);
  fIbook0 = 0;
  fNused.store(0, std::memory_order_relaxed);
  fDispatched.clear(std::memory_order_release);
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
void GeantBasket::Recycle(GeantTaskData *td) {
  // Recycle the basket to the volume scheduler.
  fManager->RecycleBasket(this, td);
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
    : TGeoExtension(), fScheduler(sch), fVolume(vol), fNumber(number), fBcap(0), fQcap(32), fActive(kFALSE),
<<<<<<< HEAD
      fCollector(collector), fThreshold(0), fNbaskets(0), fNused(0), fIbook(0), fCBasket(0), fFeeder(0),
      fDispatchList() {
=======
      fCollector(collector), fThreshold(0), fNbaskets(0), fNused(0), fCBasket(0), fLock(), fQLock(), fFeeder(0),
      fMutex() {
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  // Constructor
  fBcap = GeantPropagator::Instance()->fMaxPerBasket + 1;
  // The line below to be removed when the automatic activation schema in place
  if (collector)
    Activate();
}

//______________________________________________________________________________
GeantBasketMgr::~GeantBasketMgr() {
  // Clean up
  delete GetCBasket();
}

//______________________________________________________________________________
void GeantBasketMgr::Activate() {
  // Activate the manager for generating baskets.
  if (fActive)
    return;
  GeantBasket *basket;
  basket = new GeantBasket(fBcap, this);
  SetCBasket(basket);
  if (fCollector) {
    basket->SetMixed(true);
<<<<<<< HEAD
    //    Printf("Created collector basket manager");
=======
    Printf("Created collector basket manager");
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  }
  fActive = true;
}

//______________________________________________________________________________
<<<<<<< HEAD
GeantBasket *GeantBasketMgr::BookBasket(GeantTaskData *td) {
  // The method books a basket for track addition. It checks that
  // the basket was not replaced by other thread just after booking, in which
  // case it retries. Booking of a valid slot is eventually guaranteed.

  // Do a read from current
  GeantBasket *basket = 0;
  Int_t nbooked;
  while (1) {
    // Read current basket and book it
    basket = fCBasket.load(std::memory_order_relaxed);
    basket->fNused.fetch_add(1, std::memory_order_relaxed);
    // The copying counter is incremented before starting copying, decremented after
    basket->fNcopying.fetch_add(1, std::memory_order_relaxed);
    // The booking counter can only increment
    size_t ibook = fIbook.load(std::memory_order_relaxed);
    // Make sure baset matches ibook (not replaced)
    if (basket->fIbook0 != ibook) {
      // Basket not available anymore
      basket->fNcopying.fetch_sub(1, std::memory_order_relaxed);
      basket->fNused.fetch_sub(1, std::memory_order_relaxed);
      continue;
    }
    // Validate booking only if slot is within threshold
    nbooked = basket->fNbooked.fetch_add(1, std::memory_order_relaxed) + 1;
    if (nbooked <= basket->fThreshold)
      return basket;
    // At this point booked slot is out of range
    // No slot available, release copying
    basket->fNcopying.fetch_sub(1, std::memory_order_relaxed);
    // Any thread attempts to replace the basket corresponding to ibook here
    ReplaceBasketWeak(ibook, td);
    // Release completely basket
    basket->fNused.fetch_sub(1, std::memory_order_relaxed);
=======
GeantBasket *GeantBasketMgr::StealAndReplace(atomic_basket &current, GeantTaskData *td) {
  // Steal the current pointer content and replace with a new basket from the
  // pool. If the operation succeeds, returns the released basket which is now
  // thread local
  GeantBasket *newb, *oldb;
  // Backup cbasket. If a steal happened it does not matter
  newb = GetNextBasket(td); // backup new basket to avoid leak
  oldb = current.load();    // both local vars, so identical content
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
    newb->Recycle(td);
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  }
}

//______________________________________________________________________________
<<<<<<< HEAD
bool GeantBasketMgr::ReplaceBasketWeak(size_t expected, GeantTaskData *td) {
  // Try to replace the basket if the index is matching the expected value.
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  if (basket->fIbook0 != expected)
    return false;
  // The basket is matching the expected index,
  // now try to replace with new basket
  GeantBasket *newb = GetNextBasket(td);
  newb->fIbook0 = expected + basket->GetThreshold();
  bool replaced_by_me = fIbook.compare_exchange_weak(expected, newb->fIbook0, std::memory_order_relaxed);
  if (replaced_by_me) {
    fCBasket.store(newb, std::memory_order_relaxed);
  } else {
    td->RecycleBasket(newb);
=======
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
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  }
  return replaced_by_me;
}

//______________________________________________________________________________
<<<<<<< HEAD
bool GeantBasketMgr::ReplaceBasketStrong(size_t expected, GeantTaskData *td) {
  // Try to replace the basket if the index is matching the expected value.
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  if (basket->fIbook0 != expected)
    return false;
  // The basket is matching the expected index,
  // now try to replace with new basket
=======
Bool_t GeantBasketMgr::StealMatching(atomic_basket &global, GeantBasket *content, GeantTaskData *td) {
  // Steal the global basket if it has the right matching content
  // Prepare replacement
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  GeantBasket *newb = GetNextBasket(td);
  newb->fIbook0 = expected + basket->GetThreshold();
  bool replaced_by_me = fIbook.compare_exchange_strong(expected, newb->fIbook0);
  if (replaced_by_me) {
    fCBasket.store(newb, std::memory_order_relaxed);
  } else {
    td->RecycleBasket(newb);
  }
  return replaced_by_me;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  Int_t ncopied = 0;
  GeantBasket *basket = 0;

  // 1. Booking phase. Increments copying and booking counters. A valid basket returned
  basket = BookBasket(td);

  // 2. Track adding phase - basket not pushable
  basket->AddTrack(trackv, itr);

  // 3. Release copying so that basket can be pushed by garbage collector if need be
  basket->fNcopying.fetch_sub(1, std::memory_order_relaxed);
  // Update ncopied counter
  ncopied = basket->fNcopied.fetch_add(1, std::memory_order_relaxed) + 1;

  // 4. Last thread that completed the copy AND all slots are booked
  // has to submit the basket. The condition should ALWAYS be true for a single thread
  if (ncopied >= basket->GetThreshold()) {
    ReplaceBasketStrong(basket->fIbook0, td);
    // We can safely push the filled basket
    if (basket->TryDispatch()) {
      // Release basket
      basket->fNused.fetch_sub(1, std::memory_order_relaxed);
      Push(basket, priority, td);
      return 1;
    }
  }
  // Release basket
  basket->fNused.fetch_sub(1, std::memory_order_relaxed);
  return 0;
}

//______________________________________________________________________________
<<<<<<< HEAD
=======
Int_t GeantBasketMgr::AddTrackSingleThread(GeantTrack_v &trackv, Int_t itr, Bool_t priority, GeantTaskData *td) {
  // Copy directly from a track_v a track to the basket manager. It is
  // assumed that this manager is only handled by a single thread.
  GeantBasket *cbasket = GetCBasket();
  cbasket->GetInputTracks().AddTrack(trackv, itr);
  if (cbasket->GetNinput() >= cbasket->GetThreshold()) {
    fFeeder->push(cbasket, priority);
    SetCBasket(GetNextBasket(td));
    return 1;
  }
  return 0;
}

//______________________________________________________________________________
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
Int_t GeantBasketMgr::AddTrack(GeantTrack &track, Bool_t priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  Int_t ncopied = 0;
  GeantBasket *basket = 0;

  // 1. Booking phase. Increments copying and booking counters. A valid basket returned
  basket = BookBasket(td);

  // 2. Track adding phase - basket not pushable
  basket->AddTrack(track);

  // 3. Release copying so that basket can be pushed by garbage collector if need be
  basket->fNcopying.fetch_sub(1, std::memory_order_relaxed);
  // Update ncopied counter
  ncopied = basket->fNcopied.fetch_add(1, std::memory_order_relaxed) + 1;

  // 4. Last thread that completed the copy AND all slots are booked
  // has to submit the basket. The condition should ALWAYS be true for a single thread
  if (ncopied >= basket->GetThreshold()) {
    ReplaceBasketStrong(basket->fIbook0, td);
    // We can safely push the filled basket
    if (basket->TryDispatch()) {
      // Release basket
      basket->fNused.fetch_sub(1, std::memory_order_relaxed);
      Push(basket, priority, td);
      return 1;
<<<<<<< HEAD
    }
  }
  // Release basket
  basket->fNused.fetch_sub(1, std::memory_order_relaxed);
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrackSingleThread(GeantTrack_v &trackv, Int_t itr, Bool_t priority, GeantTaskData *td) {
  // Copy directly from a track_v a track to the basket manager. It is
  // assumed that this manager is only handled by a single thread.
  GeantBasket *cbasket = GetCBasket();
  cbasket->GetInputTracks().AddTrack(trackv, itr);
  if (cbasket->GetNinput() >= cbasket->GetThreshold()) {
    assert(cbasket->TryDispatch());
    Push(cbasket, priority, td);
    SetCBasket(GetNextBasket(td));
    return 1;
=======
    } else
      return 0;
  }
  oldb->UnLockAddingOp();
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracksNew(GeantBasketMgr *gc, GeantTaskData *td) {
  // Garbage collect tracks from the given event range. The basket gc should
  // be thread local
  if (!GetCBasket()->GetNinput())
    return 0;
  GeantBasket *cbasket = 0;
  GeantPropagator *propagator = GeantPropagator::Instance();
  while (cbasket == 0)
    cbasket = StealAndReplace(fCBasket, td);
  // cbasket is now thread local
  // Loop all tracks in the current basket and mark for copy the ones belonging
  // to the desired events
  GeantTrack_v &tracks = cbasket->GetInputTracks();
  Int_t ntracks = tracks.GetNtracks();
  Int_t icoll = 0;
  for (Int_t itr = 0; itr < ntracks; itr++) {
    if (propagator->fEvents[tracks.fEvslotV[itr]]->IsPrioritized()) {
      tracks.MarkRemoved(itr);
      icoll++;
    }
  }
  if (icoll)
    tracks.Compact(&gc->GetCBasket()->GetInputTracks());
  return icoll;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::CollectPrioritizedTracks(Int_t evmin, Int_t evmax, GeantTaskData *td) {
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
    cbasket = StealAndReplace(fCBasket, td);
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
    cbasket->Recycle(td);
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
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  }
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::GarbageCollect(GeantTaskData *td) {
  // Copy all priority tracks to the current basket and flush to queue
  // We want to steal fCBasket
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  basket->fNused.fetch_add(1, std::memory_order_relaxed);
  Int_t nbooked = basket->fNbooked.load(std::memory_order_relaxed);
  if (nbooked) {
    ReplaceBasketStrong(basket->fIbook0, td);
    if (basket->TryDispatch()) {
      // Spin till copying finished
      while (basket->fNcopying.load(std::memory_order_relaxed))
        ;
      // Release basket
      basket->fNused.fetch_sub(1, std::memory_order_relaxed);
      Push(basket, false, td);
      return 1;
    }
  }
  // Release basket
  basket->fNused.fetch_sub(1, std::memory_order_relaxed);
  return 0;
}

//______________________________________________________________________________
void GeantBasketMgr::CreateEmptyBaskets(Int_t nbaskets, GeantTaskData *td) {
  // Creates new basket for this manager
  for (auto i = 0; i < nbaskets; ++i) {
    GeantBasket *next = new GeantBasket(fBcap, this);
    if (fCollector)
      next->SetMixed(kTRUE);
    fNbaskets++;
    next->SetThreshold(fThreshold.load(std::memory_order_relaxed));
    td->RecycleBasket(next);
  }
}

//______________________________________________________________________________
void GeantBasketMgr::Push(GeantBasket *basket, Bool_t priority, GeantTaskData *td) {
  // Called whenever a basket has to be pushed to the queue. Recalculates
  // threshold for the basket manager.
  const Int_t nthreads = td->fNthreads;
  Int_t threshold = fThreshold.load(std::memory_order_relaxed);
  Int_t threshold_new = threshold * fNused.load(std::memory_order_relaxed) / nthreads;
  if ((!fCollector) && (threshold_new < fBcap) && ((threshold_new - threshold) > 4)) {
    Int_t remainder = threshold_new % 4;
    if (remainder > 0)
      threshold_new += 4 - remainder;
    fThreshold.store(threshold_new, std::memory_order_relaxed);
  }
  fNused++;
  fFeeder->push(basket, priority);
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket(GeantTaskData *td) {
  // Returns next empy basket if any available, else create a new basket.
  GeantBasket *next = td->GetNextBasket();
  if (!next) {
    next = new GeantBasket(fBcap, this);
<<<<<<< HEAD
    fNbaskets++;
=======
    if (fCollector)
      next->SetMixed(kTRUE);
    fNbaskets++;
    fNused++;
  } else {
    if (fCollector)
      next->SetMixed(kTRUE);
    else
      next->SetBasketMgr(this);
    fNused++;
>>>>>>> GEANT-133 Replacement of ROOT Materials completed -- but it just compiles
  }
  if (fCollector)
    next->SetMixed(kTRUE);
  else
    next->SetBasketMgr(this);
  next->SetThreshold(fThreshold.load(std::memory_order_relaxed));
  return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b, GeantTaskData *td) {
  // Recycle a basket.
  b->Clear();
  Int_t nthreads = td->fNthreads;
  Int_t threshold = fThreshold.load(std::memory_order_relaxed);
  Int_t threshold_new = threshold * fNused.load(std::memory_order_relaxed) / nthreads;
  if ((!fCollector) && (threshold_new > 4) && ((threshold_new - threshold) < -4)) {
    Int_t remainder = threshold_new % 4;
    threshold_new -= remainder;
    fThreshold.store(threshold_new, std::memory_order_relaxed);
  }
  td->RecycleBasket(b);
  fNused--;
}
//______________________________________________________________________________
void GeantBasketMgr::CleanBaskets(Int_t ntoclean, GeantTaskData *td) {
  // Clean a number of recycled baskets to free some memory
  Int_t ncleaned = td->CleanBaskets(ntoclean);
  fNbaskets -= ncleaned;
}

//______________________________________________________________________________
void GeantBasketMgr::Print(Option_t *) const {
  // Print info about the basket content.
  Printf("Bsk_mgr %s: current: in=%d out=%d", GetName(), GetCBasket()->GetNinput(), GetCBasket()->GetNoutput());
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
