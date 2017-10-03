#include "GeantEventServer.h"

#include "globals.h"
#include "Geant/Error.h"

#include "VBconnector.h"
#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantRunManager.h"
#include "LocalityManager.h"
#include "PrimaryGenerator.h"
#include "GeantTaskData.h"
#include "GeantBasket.h"
#include "Basket.h"
#include "StackLikeBuffer.h"
#include "MCTruthMgr.h"

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/SimpleNavigator.h"
#include "volumes/PlacedVolume.h"
#include "management/GeoManager.h"
#else
#ifdef USE_ROOT
#include "TGeoVolume.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#endif
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

using namespace vecgeom;

//______________________________________________________________________________
GeantEventServer::GeantEventServer(int nactive_max, GeantRunManager *runmgr)
  :fNevents(0), fNactiveMax(nactive_max), fNactive(0), fNserved(0), fLastActive(-1), fCurrentEvent(0),
   fNload(0), fNstored(0), fNcompleted(0), fRunMgr(runmgr),
   fFreeSlots(AdjustSize(runmgr->GetConfig()->fNbuff)),
   fPendingEvents(4096)
{
// Constructor
  assert(nactive_max > 0 && nactive_max < 4096);
  fLastActive.store(-1);
  for (size_t slot = 0; slot < size_t(fNactiveMax); ++slot)
    fFreeSlots.enqueue(slot);
  
  if (fRunMgr->GetConfig()->fRunMode == GeantConfig::kGenerator)
    fNevents = fRunMgr->GetConfig()->fNtotal;
  
  fEvents = new GeantEvent*[fNactiveMax];
  for (int i=0; i<fNactiveMax; ++i)
    fEvents[i] = nullptr;
  fGenLock.clear(std::memory_order_release);
}

//______________________________________________________________________________
GeantEventServer::~GeantEventServer()
{
// Destructor
  for (int i = 0; i < fNactiveMax; ++i) {
    delete fEvents[i];
  }
  delete [] fEvents;
}

//______________________________________________________________________________
GeantEvent *GeantEventServer::GenerateNewEvent(GeantEvent *event, GeantTaskData *td)
{
// Generates a new event in standalone GeantV mode.
  if (!fRunMgr->GetPrimaryGenerator())
    return nullptr;
  // The method has to be locked since thread safety not yet required in the generator.
  // Policy: if someone else working here, just return
  if (fGenLock.test_and_set(std::memory_order_acquire)) return nullptr;
  // Now just get next event from the generator
  GeantEventInfo eventinfo = fRunMgr->GetPrimaryGenerator()->NextEvent();
  int ntracks = eventinfo.ntracks;
  if (!ntracks) return nullptr;
  
  if (!event)
    event = new GeantEvent();
  else
    event->Clear();
  fNload++;
  event->SetNprimaries(ntracks);
  event->SetVertex(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);

  // Initialize navigation path for the vertex
  Volume_t *vol = 0;
  // Initialize the start path
  VolumePath_t *startpath = VolumePath_t::MakeInstance(fRunMgr->GetConfig()->fMaxDepth);
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::SimpleNavigator nav;
  startpath->Clear();
  nav.LocatePoint(GeoManager::Instance().GetWorld(),
                  Vector3D<Precision>(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert), *startpath, true);
  vol = const_cast<Volume_t *>(startpath->Top()->GetLogicalVolume());
  VBconnector *link = static_cast<VBconnector *>(vol->GetBasketManagerPtr());
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
  TGeoNode *geonode = nav->FindNode(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);
  vol = geonode->GetVolume();
  VBconnector *link = static_cast<VBconnector *>(vol->GetFWExtension());
  startpath->InitFromNavigator(nav);
#endif
  // There are needed for v2 only
  fBindex = link->index;
  td->fVolume = vol;
  
  // Populate event with primary tracks
  for (int itr=0; itr<ntracks; ++itr) {
    GeantTrack &track = td->GetNewTrack();
    track.fParticle = event->AddPrimary(&track);
    track.SetPrimaryParticleIndex(itr); 
    track.SetPath(startpath);
    track.SetNextPath(startpath);
    fRunMgr->GetPrimaryGenerator()->GetTrack(itr, track);
    if (!track.IsNormalized())
      track.Print("Not normalized");
    track.fBoundary = false;
    track.fStatus = kNew;
    event->fNfilled++;
    if (fRunMgr->GetMCTruthMgr()) fRunMgr->GetMCTruthMgr()->AddTrack(track);
  }
 
  fGenLock.clear(std::memory_order_release);
  if (!AddEvent(event)) return nullptr;
  return event;
}

//______________________________________________________________________________
bool GeantEventServer::AddEvent(GeantEvent *event)
{
// Adds one event into the queue of pending events.
  int evt = fNload.fetch_add(1);

  if (fRunMgr->GetConfig()->fRunMode == GeantConfig::kExternalLoop) {
    // The vertex must be defined
    vecgeom::Vector3D<double> vertex = event->GetVertex();
    int ntracks = event->GetNprimaries();
    
    // Initialize navigation path for the vertex
    Volume_t *vol = 0;
    // Initialize the start path
    VolumePath_t *startpath = VolumePath_t::MakeInstance(fRunMgr->GetConfig()->fMaxDepth);
#ifdef USE_VECGEOM_NAVIGATOR
    vecgeom::SimpleNavigator nav;
    startpath->Clear();
    nav.LocatePoint(GeoManager::Instance().GetWorld(), vertex, *startpath, true);
    vol = const_cast<Volume_t *>(startpath->Top()->GetLogicalVolume());
    VBconnector *link = static_cast<VBconnector *>(vol->GetBasketManagerPtr());
#else
    TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
    if (!nav)
      nav = gGeoManager->AddNavigator();
    TGeoNode *geonode = nav->FindNode(vertex.x(), vertex.y(), vertex.z());
    vol = geonode->GetVolume();
    VBconnector *link = static_cast<VBconnector *>(vol->GetFWExtension());
    startpath->InitFromNavigator(nav);
#endif
    // There are needed for v2 only
    fBindex = link->index;
    //td->fVolume = vol;
    // Check and fix tracks
    for (int itr=0; itr<ntracks; ++itr) {
      GeantTrack &track = *event->GetPrimary(itr);
      track.SetPrimaryParticleIndex(itr); 
      track.SetPath(startpath);
      track.SetNextPath(startpath);
      track.SetEvent(evt);
      if (!track.IsNormalized())
        track.Print("Not normalized");
      track.fBoundary = false;
      track.fStatus = kNew;
      event->fNfilled++;
      if (fRunMgr->GetMCTruthMgr()) fRunMgr->GetMCTruthMgr()->AddTrack(track);
    }
    // Release path object
    VolumePath_t::ReleaseInstance(startpath);
  }
    
  if (!fPendingEvents.enqueue(event)) {
    Error("AddEvent", "Event pool is full");
    return false;
  }
  event->SetPriorityThr(fRunMgr->GetConfig()->fPriorityThr);
  // Update number of stored events
  int nstored = fNstored.fetch_add(1) + 1;
  if (nstored <= fNactiveMax) fNtracksInit += event->GetNprimaries();
  if (nstored == fNactiveMax) {
    // Check consistency of parallelism settings
    int nthreads = fRunMgr->GetNthreadsTotal();
    int basket_size = fRunMgr->GetConfig()->fNperBasket;
    int ntracksperevent = fNtracksInit/fNactiveMax;
    fNbasketsInit = fNtracksInit/basket_size;
    if (fNtracksInit % basket_size > 0) fNbasketsInit++;
    printf("=== Imported %d primaries from %d buffered events\n", fNtracksInit, fNactiveMax);
    printf("=== Buffering %d baskets of size %d feeding %d threads\n", fNbasketsInit, basket_size, nthreads);
    if (fNbasketsInit < nthreads) {
      printf("### WARNING!    Concurrency settings are not optimal. Not enough baskets to feed all threads.\n###\n");
      printf("###                          Rule of thumb:\n");
      printf("###  ==================================================================\n");
      printf("### ||  Nbuff_events * Nprimaries_per_event > Nthreads * basket_size  ||\n");
      printf("###  ==================================================================\n###\n");
      printf("### Either of the following options can solve this problem:\n");
      printf("###    => increase number of buffered events slots to minimum %d\n", nthreads * basket_size / ntracksperevent);
      printf("###    => decrease number of threads to maximum %d\n", vecCore::math::Max(1, fNtracksInit/basket_size));    
    }
  }

  fEventsServed = false;
  return true;
}

//______________________________________________________________________________
GeantTrack *GeantEventServer::GetNextTrack(unsigned int &error)
{
// Fetch next track of the current event. Increments current event if no more
// tracks. If current event matches last activated one, resets fHasTracks flag.
// If max event fully dispatched, sets the fDone flag.

  GeantEvent *event = nullptr;
  int itr;
  while (1) {
    // Book next track index
    bool valid;
    event = fEvent.load();
    itr = event->DispatchTrack(valid);
    if (!valid) {
      assert(event->IsDispatched());
      // Current event dispatched, try to fetch from next event
      event = ActivateEvent(event, error);
      if (!event) {
        if (error) return nullptr;
        continue;
      }
      itr = event->DispatchTrack(valid);
    }
    if (valid) break;
  }
  GeantTrack *track = event->GetPrimary(itr);
  track->fEvent = event->GetEvent();
  track->fEvslot = event->GetSlot();
  return track;
}

//______________________________________________________________________________
int GeantEventServer::FillBasket(GeantTrack_v &tracks, int ntracks, unsigned int &error)
{
// Fill concurrently a basket of tracks, up to the requested number of tracks.
// The error codes returned: 0-normal fill 1-partial fill 2-
// The client should test first the track availability using HasTracks().
  if (!fHasTracks) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack(error);
    if (!track) break;
    tracks.AddTrack(*track);
    ndispatched++;
  }
  if (fInitialPhase) {
    int nserved = fNserved.fetch_add(1) + 1;
    if (nserved >= fNbasketsInit) fInitialPhase = false;
  }
  return ndispatched;
}

//______________________________________________________________________________
int GeantEventServer::FillBasket(Basket *basket, int ntracks, unsigned int &error)
{
// Fill concurrently a basket of tracks, up to the requested number of tracks.
// The client should test first the track availability using HasTracks().
  if (!fHasTracks) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack(error);
    if (!track) break;
    basket->AddTrack(track);
    ndispatched++;
  }
  if (fInitialPhase) {
    int nserved = fNserved.fetch_add(1) + 1;
    if (nserved >= fNbasketsInit) fInitialPhase = false;
  }
  return ndispatched;
}

//______________________________________________________________________________
int GeantEventServer::FillStackBuffer(StackLikeBuffer *buffer, int ntracks, unsigned int &error)
{
// Fill concurrently up to the requested number of tracks into a stack-like buffer.
// The client should test first the track availability using HasTracks().

// *** I should template on the container to be filled, making sure that all
//     containers provide AddTrack(GeantTrack *)
  if (!fHasTracks) {
    error = kDone;
    return 0;
  }
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack(error);
    if (!track) break;
    buffer->AddTrack(track);
    ndispatched++;
  }
  if (fInitialPhase) {
    int nserved = fNserved.fetch_add(1) + 1;
    if (nserved >= fNbasketsInit) fInitialPhase = false;
  }
  
  if (ndispatched > 0) error = kNoerr;
  return ndispatched;
}

//______________________________________________________________________________
GeantEvent *GeantEventServer::ActivateEvent(GeantEvent *event, unsigned int &error)
{
// Activates one event replacing the current one (if matching the expected value).
// The method can fail due to one of the following reasons:
// - All events already served (nactive = ntotal) in generator mode
// - All events from buffer served (nactive = nstored) in external loop mode
// - No slots available (should never happen since activation is driven by slot release)
// - All buffered events were served (should not happen in generator mode)
  if (fEventsServed) {
    error = kDone;
    return nullptr;
  }
  size_t slot = 0;

  // A slot should always be available when calling ActivateEvent
  if (!fFreeSlots.dequeue(slot)) {
    error = kCSlots;
    return nullptr;
  }
  printf("(A) slot #%ld taken\n", slot);

  // The event stored at the retrieved slot has to be transported
  if (fEvents[slot]) { assert(fEvents[slot]->Transported()); }
  
  // An event should also be prefetched into the waiting queue. If none
  // is available it is due to generator being too slow.
  if (!fPendingEvents.dequeue(fEvent, event)) {
    fFreeSlots.enqueue(slot);
    if (!fPendingEvents.size()) error |= kCEvents;
    else error = kNoerr;
    printf("(A) slot #%ld released (error %d)\n", slot, error);
    return nullptr;
  }
  // Activate the event
  fEvents[slot] = event;
  fEvents[slot]->SetSlot(slot);
  int nactive = fNactive.fetch_add(1) + 1;
  fEvents[slot]->SetEvent(nactive - 1);
  // Check if all events were served
  if (fRunMgr->GetConfig()->fRunMode == GeantConfig::kGenerator) {
    if (nactive == fNevents) fEventsServed = true;
  } else {
    if (nactive == fNstored.load()) fEventsServed = true;
  }
  fHasTracks = true;
  error = kNoerr;
  // At this point other threads can feed from this new event.
  // printf("Activated event %p on slot %ld\n", event, slot);
  return event;
}

//______________________________________________________________________________
void GeantEventServer::CompletedEvent(GeantEvent *event, GeantTaskData *td)
{
// Signals that event 'evt' was fully transported.
  size_t slot = event->GetSlot();
  fNcompleted++;
  assert(event->Transported());
  assert(event == fEvents[slot]);
  fEvents[slot] = nullptr;
  if (fRunMgr->GetConfig()->fRunMode == GeantConfig::kGenerator &&
      fNstored.load() < fNevents)
    GenerateNewEvent(event, td);
  // Now we can release the slot, making sure the current event has changed
  while (!fEventsServed && event == fEvent.load()) {}
  fFreeSlots.enqueue(slot);
  printf("(C) slot #%ld released by task %d\n", slot, td->fTid);
}

} // GEANT_IMPL_NAMESPACE
} // Geant
