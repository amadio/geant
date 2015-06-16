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
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fNstalled(0), fLock(), 
      fDispatched(), fReplaced(false), fThreshold(0), fTracksIn(), fTracksOut() {
  // Dummy constructor.
  fLock.clear();
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, GeantBasketMgr *mgr)
    : TObject(), fManager(mgr), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fNstalled(0), fLock(), fDispatched(), 
      fReplaced(false), fThreshold(size),
      fTracksIn(size, GeantPropagator::Instance()->fMaxDepth),
      fTracksOut(size, GeantPropagator::Instance()->fMaxDepth) {
  // Default constructor.
  if (!mgr->GetVolume() || mgr->IsCollector()) 
    SetMixed(true);
  fLock.clear();
}

//______________________________________________________________________________
GeantBasket::GeantBasket(Int_t size, Int_t depth)
    : TObject(), fManager(0), fNcopying(0), fNbooked(0), fNcopied(0), fNused(0), fNstalled(0), fLock(),
      fDispatched(), fReplaced(false), fThreshold(size), fTracksIn(size, depth), fTracksOut(size, depth) {
  // Default constructor.
  fLock.clear();
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
  lock();
  SetMixed(fManager->IsCollector());
  fTracksIn.Clear(option);
  fTracksOut.Clear(option);
  fNbooked.store(0);
  fNcopying.store(0);
  fNcopied.store(0);
  fNused.store(0);
  fReplaced.store(false);
  fDispatched.clear(std::memory_order_release);
  fReplaced.store(false);
  unlock();
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
      fActive(kFALSE), fCollector(collector), fThreshold(0), fNbaskets(0), fNused(0), 
      fCBasket(0), fDLock(), fStealLock(), fFeeder(0), fDispatchList() {
  // Constructor
  fDLock.clear();
  fStealLock.clear();
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
GeantBasket *GeantBasketMgr::BookBasket(atomic_basket &current, GeantTaskData *td,
                                        Int_t &nbooked, Int_t &ncopying) {
  // The method books a basket for track addition. It checks that
  // the basket was not replaced by other thread just after booking, in which 
  // case it retries. Booking of a valid slot is eventually guaranteed.
  
  // Do a read from current
  GeantBasket *basket = 0;
  while (1) {
    // Read current basket in non-replace mode and mark as copying (busy)
    lock_steal();
    basket = current.load();
    basket->fNused++;
    // Book a slot for copying.
    ncopying = ++basket->fNcopying;
    nbooked = ++basket->fNbooked;
    unlock_steal();
    // At this point the booking slot may be greater that the threshold
    // Validate booking only if slot is within threshold
    if (nbooked <= basket->fThreshold) return basket;
//      if (nbooked == 1) AddToDispatch(basket);
    // At this point booked slot is out of range
    // No slot available, release copying
    ncopying = --basket->fNcopying;
    // The first one who finds an invalid slot tries to replace basket in none else did it
    // The statement below should stand for at least one. Managing the
    // replacement is NOT critical here, but just an optimisation 
    if ( nbooked > basket->fThreshold )
      ReplaceBasket(current, basket, false, td);
    // Release completely basket
    basket->fNused--;
  }
}

//______________________________________________________________________________
bool GeantBasketMgr::ReplaceBasket(atomic_basket &current, GeantBasket *expected, bool strong, 
                                   GeantTaskData *td)
{
  // Replace content of current with a new basket, if not already done. Returns true only
  // if replacement was performed by this.
  if (expected->fReplaced.load() || (current.load() != expected)) return false;
  GeantBasket *newb = GetNextBasket(td);
  bool replaced_by_me = false;
  lock_steal();
  // The compare exchange below can be intruding into setting copying busy flag - prevent that
  if (strong) replaced_by_me = fCBasket.compare_exchange_strong(expected, newb);
  else        replaced_by_me = fCBasket.compare_exchange_weak(expected, newb);
  if (replaced_by_me) 
    expected->fReplaced.store(true);
  unlock_steal();
  if (!replaced_by_me) newb->Recycle(td);
  return replaced_by_me;
}   
  
//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack_v &trackv, Int_t itr, Bool_t priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  Int_t nbooked = 0;
  Int_t ncopying = 0;
  Int_t ncopied = 0;
  GeantBasket *basket = 0;

  // 1. Booking phase. Increments copying counter and number of slots. Anyone
  //    can steal the basket, but not push it.
  basket = BookBasket(fCBasket, td, nbooked, ncopying);

  // 2. Track adding phase - basket not pushable
  basket->AddTrack(trackv, itr);
//  Int_t ncopied = basket->TrackCopied();
  // Write to copying counter and reading from booking counter must be coupled 
#ifdef BASKETIZER_DEBUG
  assert(ncopying > 0);
  assert(nbooked > 0);
#endif  
  basket->fNcopying--;
  // Update ncopied counter
  ncopied = basket->fNcopied.load(std::memory_order_relaxed);
  while (!basket->fNcopied.compare_exchange_weak(ncopied, ncopied+1, std::memory_order_relaxed)) ;
  ncopied++; // now contains the new value when changed by THIS thread

  // 3. Last thread that completed the copy AND all slots are booked
  // has to submit the basket. The condition should ALWAYS be true for a single thread
  if (ncopied >= basket->GetThreshold()) {
    ReplaceBasket(fCBasket, basket, true, td);
#ifdef BASKETIZER_DEBUG
//    assert(basket->fReplaced.load());
#endif
    // We can safely push the filled basket
    if (basket->TryDispatch()) {
//      RemoveFromDispatch(basket);
//      basket->Clear();
      // Spin till copying finished (should not spin)
      while (basket->fNcopying.load()) ;
      // Release basket
      basket->fNused--;
      fFeeder->push(basket, priority);
      return 1;
    }
  }
  // Release basket
  basket->fNused--;
  return 0;
}

//______________________________________________________________________________
Int_t GeantBasketMgr::AddTrack(GeantTrack &track, Bool_t priority, GeantTaskData *td) {
  // Add a track to the current basket. If the track number reaches the
  // threshold, the basket is added to the feeder queue and replaced by an empty
  // one. Returns the number of dispatched baskets
  // Has to work concurrently
  Int_t nbooked = 0;
  Int_t ncopying = 0;
  Int_t ncopied = 0;
  GeantBasket *basket = 0;

  // 1. Booking phase. Increments copying counter and number of slots. Anyone
  //    can steal the basket, but not push it.
  basket = BookBasket(fCBasket, td, nbooked, ncopying);

  // 2. Track adding phase - basket not pushable
  basket->AddTrack(track);
//  Int_t ncopied = basket->TrackCopied();
  // Write to copying counter and reading from booking counter must be coupled 
#ifdef BASKETIZER_DEBUG
  assert(ncopying > 0);
  assert(nbooked > 0);
#endif
  basket->fNcopying--;
  // Update ncopied counter
  ncopied = basket->fNcopied.load(std::memory_order_relaxed);
  while (!basket->fNcopied.compare_exchange_weak(ncopied, ncopied+1, std::memory_order_relaxed)) ;
  ncopied++; // now contains the new value when changed by THIS thread

  // 3. Last thread that completed the copy AND all slots are booked
  // has to submit the basket. The condition should ALWAYS be true for a single thread
  if (ncopied >= basket->GetThreshold()) {
    ReplaceBasket(fCBasket, basket, true, td);
#ifdef BASKETIZER_DEBUG
//    assert(basket->fReplaced.load());
#endif
    // We can safely push the filled basket
    if (basket->TryDispatch()) {
//      RemoveFromDispatch(basket);
//      basket->Clear();
      // Spin till copying finished (should not spin)
      while (basket->fNcopying.load()) ;
      // Release basket
      basket->fNused--;
      fFeeder->push(basket, priority);
      return 1;
    }
  }
  // Release basket
  basket->fNused--;
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
    fFeeder->push(cbasket, priority);    
    SetCBasket(GetNextBasket(td));
    return 1;
  }
  return 0;
}
  
//______________________________________________________________________________
Int_t GeantBasketMgr::GarbageCollect(GeantTaskData *td) {
  // Copy all priority tracks to the current basket and flush to queue
  // We want to steal fCBasket
//  CheckStalled();
//  int nstalled = fDispatchList.size();
//  if (nstalled) Printf(" === %s: %d baskets stalled", GetName(), nstalled);
//  for (auto basket1 : fDispatchList) Printf("   %p: %d/%d", (void*)basket1, basket1->fNcopied.load(), basket1->GetThreshold());
  lock_steal();
  GeantBasket *basket = fCBasket.load();
  basket->fNused++;
  Int_t nbooked = basket->fNbooked.load();
  unlock_steal();
  if (nbooked) {
    ReplaceBasket(fCBasket, basket, true, td);
#ifdef BASKETIZER_DEBUG
//    assert(basket->fReplaced.load());
#endif
    if (basket->TryDispatch()) {
      // Spin till copying finished (should not spin)
      while (basket->fNcopying.load()) ;
      // Release basket
      basket->fNused--;
      fFeeder->push(basket, false);
      return 1;
    }
  }
  // Release basket
  basket->fNused--;
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
    next->SetThreshold(fThreshold.load());
    td->RecycleBasket(next);
  }  
}

//______________________________________________________________________________
GeantBasket *GeantBasketMgr::GetNextBasket(GeantTaskData *td) {
  // Returns next empy basket if any available, else create a new basket.
  GeantBasket *next = 0;
  const Int_t nthreads = td->fNthreads;
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
  next = td->GetNextBasket();
  if (!next) {
    next = new GeantBasket(fBcap, this);
    if (fCollector) next->SetMixed(kTRUE);
    fNbaskets++;
    fNused++;
  } else {
    if (fCollector) next->SetMixed(kTRUE);
    else            next->SetBasketMgr(this);
    fNused++;
  }
  next->SetThreshold(fThreshold.load());
  return next;
}

//______________________________________________________________________________
void GeantBasketMgr::RecycleBasket(GeantBasket *b, GeantTaskData *td) {
  // Recycle a basket.
  b->Clear();
  Int_t nthreads = td->fNthreads;
  Int_t threshold = fThreshold.load();
  Int_t threshold_new = threshold * fNused.load() / nthreads;
  if ((!fCollector) && (threshold_new > 4) && ((threshold_new - threshold) < -4)) {
    Int_t remainder = threshold_new % 4;
    threshold_new -= remainder;
    fThreshold.store(threshold_new);
    //      Printf("volume %s:      ntotal=%d   nused=%d   thres=%d    size=%ld", GetName(),
    //      fNbaskets.load(), fNused.load(), threshold_new, Sizeof());
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
