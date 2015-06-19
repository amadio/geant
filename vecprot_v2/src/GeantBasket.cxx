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
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNavigator.h"

ClassImp(GeantBasket)

//______________________________________________________________________________
GeantBasket::GeantBasket()
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), 
      fDispatched(), fThreshold(0), fTracksIn(), fTracksOut() {
  // Dummy constructor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
    : TObject(), fManager(mgr), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), 
      fDispatched(), fThreshold(size),
      fTracksIn(size, GeantPropagator::Instance()->fMaxDepth),
      fTracksOut(size, GeantPropagator::Instance()->fMaxDepth) {
  // Default constructor.
  if (!mgr->GetVolume() || mgr->IsCollector()) 
    SetMixed(true);
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, Int_t depth)
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fIbook0(0), 
      fDispatched(), fThreshold(size), 
      fTracksIn(size, depth), fTracksOut(size, depth) {
  // Default constructor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket() {
  // Destructor.
}

//______________________________________________________________________________
Int_t GeantBasket::AddTrack(GeantTrack &track) {
  // Add a new track to this basket. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(track);
  return ( fTracksIn.AddTrackSync(track) );
}

//______________________________________________________________________________
Int_t GeantBasket::AddTrack(GeantTrack_v &tracks, Int_t itr) {
  // Add track from a track_v array. Has to work concurrently.

  // Activating the line below adds non-concurrently an input track
  //   fTracksIn.AddTrack(tracks, itr, kTRUE);
  return ( fTracksIn.AddTrackSync(tracks, itr) );
}

//______________________________________________________________________________
void GeantBasket::AddTracks(GeantTrack_v &tracks, Int_t istart, Int_t iend) {
  // Add multiple tracks from a track_v array.
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
    : TGeoExtension(), fScheduler(sch), fVolume(vol), fNumber(number), fBcap(0), fQcap(32),
      fActive(kFALSE), fCollector(collector), fThreshold(0), fNbaskets(0), fNused(0), fIbook(0),
      fCBasket(0), fFeeder(0), fDispatchList() {
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
void GeantBasketMgr::Activate()
{
  // Activate the manager for generating baskets.
  if (fActive) return;
  GeantBasket *basket;
  basket = new GeantBasket(fBcap, this);
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
    if (nbooked <= basket->fThreshold) return basket;
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
bool GeantBasketMgr::ReplaceBasketWeak(size_t expected,
                                   GeantTaskData *td)
{
// Try to replace the basket if the index is matching the expected value.
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  if (basket->fIbook0 != expected) return false;
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
bool GeantBasketMgr::ReplaceBasketStrong(size_t expected,
                                   GeantTaskData *td)
{
// Try to replace the basket if the index is matching the expected value.
  GeantBasket *basket = fCBasket.load(std::memory_order_relaxed);
  if (basket->fIbook0 != expected) return false;
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
  ncopied = basket->fNcopied.fetch_add(1, std::memory_order_relaxed)+1;

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
  ncopied = basket->fNcopied.fetch_add(1, std::memory_order_relaxed)+1;

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
Int_t GeantBasketMgr::AddTrackSingleThread(GeantTrack_v &trackv, Int_t itr, Bool_t priority,
                       GeantTaskData *td) {
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
      while (basket->fNcopying.load(std::memory_order_relaxed)) ;
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
void GeantBasketMgr::CreateEmptyBaskets(Int_t nbaskets, GeantTaskData *td)
{
// Creates new basket for this manager
  for (auto i=0; i<nbaskets; ++i) {
    GeantBasket *next = new GeantBasket(fBcap, this);
    if (fCollector) next->SetMixed(kTRUE);
    fNbaskets++;
    next->SetThreshold(fThreshold.load(std::memory_order_relaxed));
    td->RecycleBasket(next);
  }  
}

//______________________________________________________________________________
void GeantBasketMgr::Push(GeantBasket *basket, Bool_t priority, GeantTaskData *td)
{
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
    fNbaskets++;
  }
  if (fCollector) next->SetMixed(kTRUE);
  else            next->SetBasketMgr(this);
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
