#include "GeantBasket.h"
#include "globals.h"
#include "GeantTrackVec.h"
#include "GeantEvent.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
//#include "PhysicsProcessOld.h"
#include "WorkloadManager.h"
#ifdef USE_VECGEOM_NAVIGATOR
#else
#include "TGeoManager.h"
#include "TGeoNavigator.h"
#endif
#include "Geant/Error.h"

//______________________________________________________________________________
GeantBasket::GeantBasket()
    : fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
      fThreshold(0), fTracksIn(), fIsMixed(false) {
  // Dummy constructor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(GeantPropagator* prop, int size, GeantBasketMgr *mgr)
   : fManager(mgr), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
     fThreshold(size), fTracksIn(size, prop->fConfig->fMaxDepth), fIsMixed(false) {
  // Default constructor.
  if (!mgr->GetVolume() || mgr->IsCollector())
    SetMixed(true);
}

//______________________________________________________________________________
GeantBasket::GeantBasket(int size, int depth)
    : fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), fDispatched(),
      fThreshold(size), fTracksIn(size, depth), fIsMixed(false) {
  // Default constructor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket() {
  // Destructor.
}

//______________________________________________________________________________
int GeantBasket::AddTrack(GeantTrack &track) {
  // Add a new track to this basket. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(track);
  return (fTracksIn.AddTrackSync(track));
}

//______________________________________________________________________________
int GeantBasket::AddTrack(GeantTrack_v &tracks, int itr) {
  // Add track from a track_v array. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(tracks, itr, true);
  return (fTracksIn.AddTrackSync(tracks, itr));
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, int istart, int iend) {
  // Add multiple tracks from a track_v array.
  fTracksIn.AddTracks(tracks, istart, iend, true);
}

//______________________________________________________________________________
void GeantBasket::Clear(const char *option) {
  // Clear basket content.
  SetMixed(fManager->IsCollector());
  fTracksIn.Clear(option);
  fNbooked.store(0, std::memory_order_relaxed);
  fNcopying.store(0, std::memory_order_relaxed);
  fNcopied.store(0, std::memory_order_relaxed);
  fIbook0 = 0;
  fNused.store(0, std::memory_order_relaxed);
  fDispatched.clear(std::memory_order_release);
}

//______________________________________________________________________________
bool GeantBasket::Contains(int evstart, int nevents) const {
  // Checks if any of the tracks in the input array belongs to the given event
  // range.
  return fTracksIn.Contains(evstart, nevents);
}

//______________________________________________________________________________
Volume_t *GeantBasket::GetVolume() const {
  // Returns volume for this basket
  return fManager->GetVolume();
}

//______________________________________________________________________________
void GeantBasket::Print(const char *) const {
  // Print basket content.
  Geant::Printf("*** basket : ntracks=");
}

//______________________________________________________________________________
void GeantBasket::PrintTrack(int /*itr*/, bool /*input*/) const {
  // Print a given track.
}

//______________________________________________________________________________
void GeantBasket::Recycle(GeantTaskData *td) {
  // Recycle the basket to the volume scheduler.
  fManager->RecycleBasket(this, td);
}

//______________________________________________________________________________
void GeantBasket::SetThreshold(int threshold) {
  // Set transport threshold for the basket
  if (threshold > fThreshold) {
    if (fTracksIn.Capacity() < threshold)
      fTracksIn.Resize(threshold);
  }
  fThreshold = threshold;
}

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________
//______________________________________________________________________________
GeantBasketMgr::GeantBasketMgr(GeantPropagator* prop, GeantScheduler *sch, Volume_t *vol, int number, bool collector)
#ifdef USE_ROOT
  : TGeoExtension(), fScheduler(sch), fVolume(vol), fNumber(number), fBcap(0), fQcap(32), fActive(false),
#else
  : fScheduler(sch), fVolume(vol), fNumber(number), fBcap(0), fQcap(32), fActive(false),
#endif
    fCollector(collector), fThreshold(prop->fConfig->fNperBasket), fNbaskets(0), fNused(0), fIbook(0), fCBasket(0), fFeeder(0),
    fDispatchList() {
  // Constructor
  fBcap = prop->fConfig->fMaxPerBasket + 1;
  // The line below to be removed when the automatic activation schema in place
  if (collector)
    Activate(prop);
}

//______________________________________________________________________________
GeantBasketMgr::~GeantBasketMgr() {
  // Clean up
  delete GetCBasket();
}

//______________________________________________________________________________
void GeantBasketMgr::Activate(GeantPropagator* prop) {
  // Activate the manager for generating baskets.
  if (fActive)
    return;
  GeantBasket *basket;
  basket = new GeantBasket(prop, fBcap, this);
  basket->SetThreshold(fThreshold.load(std::memory_order_relaxed));
  SetCBasket(basket);
  if (fCollector) {
    basket->SetMixed(true);
    //    Printf("Created collector basket manager");
  }
  fActive = true;
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::BookBasket(GeantTaskData *td) {
  // The method books a basket for track addition. It checks that
  // the basket was not replaced by other thread just after booking, in which
  // case it retries. Booking of a valid slot is eventually guaranteed.

  // Do a read from current
  GeantBasket *basket = 0;
  int nbooked;
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
  }
}

//______________________________________________________________________________
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
  }
  return replaced_by_me;
}

//______________________________________________________________________________
bool GeantBasketMgr::ReplaceBasketStrong(size_t expected, GeantTaskData *td) {
  // Try to replace the basket if the index is matching the expected value.
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  if (basket->fIbook0 != expected)
    return false;
  // The basket is matching the expected index,
  // now try to replace with new basket
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
int GeantBasketMgr::AddTrack(GeantTrack_v &trackv, int itr, bool priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  int ncopied = 0;
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
int GeantBasketMgr::AddTrack(GeantTrack &track, bool priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  int ncopied = 0;
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
    }
  }
  // Release basket
  basket->fNused.fetch_sub(1, std::memory_order_relaxed);
  return 0;
}

//______________________________________________________________________________
int GeantBasketMgr::AddTrackSingleThread(GeantTrack_v &trackv, int itr, bool priority, GeantTaskData *td) {
  // Copy directly from a track_v a track to the basket manager. It is
  // assumed that this manager is only handled by a single thread.
  GeantBasket *cbasket = GetCBasket();
  cbasket->GetInputTracks().AddTrack(trackv, itr);
  if (cbasket->GetNinput() >= cbasket->GetThreshold()) {
    assert(cbasket->TryDispatch());
    Push(cbasket, priority, td);
    SetCBasket(GetNextBasket(td));
    return 1;
  }
  return 0;
}

//______________________________________________________________________________
int GeantBasketMgr::GarbageCollect(GeantTaskData *td) {
  // Copy all priority tracks to the current basket and flush to queue
  // We want to steal fCBasket
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  basket->fNused.fetch_add(1, std::memory_order_relaxed);
  int nbooked = basket->fNbooked.load(std::memory_order_relaxed);
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
void GeantBasketMgr::CreateEmptyBaskets(int nbaskets, GeantTaskData *td) {
  // Creates new basket for this manager
  for (auto i = 0; i < nbaskets; ++i) {
    GeantBasket *next = new GeantBasket(td->fPropagator,fBcap, this);
    if (fCollector)
      next->SetMixed(true);
    fNbaskets++;
    next->SetThreshold(fThreshold.load(std::memory_order_relaxed));
    td->RecycleBasket(next);
  }
}

//______________________________________________________________________________
void GeantBasketMgr::Push(GeantBasket *basket, bool priority, GeantTaskData *td) {
  // Called whenever a basket has to be pushed to the queue. Recalculates
  // threshold for the basket manager.
  const int nthreads = td->fNthreads;
  int threshold = fThreshold.load(std::memory_order_relaxed);
  int threshold_new = threshold * fNused.load(std::memory_order_relaxed) / nthreads;
  if ((!fCollector) && (threshold_new < fBcap) && ((threshold_new - threshold) > 4)) {
    int remainder = threshold_new % 4;
    if (remainder > 0)
      threshold_new += 4 - remainder;
    fThreshold.store(threshold_new, std::memory_order_relaxed);
  }
  fNused++;
  fFeeder->push_force(basket, priority);
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket(GeantTaskData *td) {
  // Returns next empy basket if any available, else create a new basket.
  GeantBasket *next = td->GetNextBasket();
  if (!next) {
    next = new GeantBasket(td->fPropagator,fBcap, this);
    fNbaskets++;
  }
  if (fCollector)
    next->SetMixed(true);
  else
    next->SetBasketMgr(this);
  next->SetThreshold(fThreshold.load(std::memory_order_relaxed));
  return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b, GeantTaskData *td) {
  // Recycle a basket.
  b->Clear();
  int nthreads = td->fNthreads;
  int threshold = fThreshold.load(std::memory_order_relaxed);
  int threshold_new = threshold * fNused.load(std::memory_order_relaxed) / nthreads;
  if ((!fCollector) && (threshold_new > 4) && ((threshold_new - threshold) < -4)) {
    int remainder = threshold_new % 4;
    threshold_new -= remainder;
    fThreshold.store(threshold_new, std::memory_order_relaxed);
  }
  td->RecycleBasket(b);
  fNused--;
}
//______________________________________________________________________________
void GeantBasketMgr::CleanBaskets(int ntoclean, GeantTaskData *td) {
  // Clean a number of recycled baskets to free some memory
  int ncleaned = td->CleanBaskets(ntoclean);
  fNbaskets -= ncleaned;
}

//______________________________________________________________________________
void GeantBasketMgr::Print(const char *) const {
  // Print info about the basket content.
  Geant::Printf("Bsk_mgr %s: current: tracks=%d", GetName(), GetCBasket()->GetNinput());
}

//______________________________________________________________________________
void GeantBasketMgr::PrintSize() const {
  // Print detailed info about size.
  size_t size = Sizeof();
  size_t sizeb = 0;
  if (GetCBasket())
    sizeb = GetCBasket()->Sizeof();
  Geant::Printf("Bsk_mgr %s: %d baskets of size %ld:    %ld bytes", GetName(), GetNbaskets(), sizeb, size);
}
