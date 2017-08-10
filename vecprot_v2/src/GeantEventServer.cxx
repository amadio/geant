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
GeantEventServer::GeantEventServer(int event_capacity, GeantRunManager *runmgr)
  :fNevents(event_capacity), fNactive(0), fNserved(0), fLastActive(-1), fCurrentEvent(0),
   fNload(0), fNstored(0), fNcompleted(0), fRunMgr(runmgr),
   fFreeSlots(AdjustSize(runmgr->GetConfig()->fNbuff)),
   fPendingEvents(4096)
{
// Constructor
  assert(event_capacity > 0);
  fLastActive.store(-1);
  fNactiveMax = runmgr->GetConfig()->fNbuff;
  for (size_t slot = 0; slot < size_t(fNactiveMax); ++slot)
    fFreeSlots.enqueue(slot);
  fEvents = new GeantEvent*[event_capacity];
  for (int i=0; i<event_capacity; ++i) {
    fEvents[i] = new GeantEvent();
    fEvents[i]->SetPriorityThr(runmgr->GetConfig()->fPriorityThr);
  }
  fGenLock.clear(std::memory_order_release);
}

//______________________________________________________________________________
GeantEventServer::~GeantEventServer()
{
// Destructor
  for (int i = 0; i < fNload.load(); ++i) {
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
  int evt = fNload.fetch_add(1);
  event->SetEvent(evt);
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
    track.SetEvent(evt);
    fRunMgr->GetPrimaryGenerator()->GetTrack(itr, track);
    if (!track.IsNormalized())
      track.Print("Not normalized");
    track.fBoundary = false;
    track.fStatus = kNew;
    event->fNfilled++;
    if (fRunMgr->GetMCTruthMgr()) fRunMgr->GetMCTruthMgr()->AddTrack(track);
  }
 
  fGenLock.clear(std::memory_order_release);
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
  // Update number of stored events
  fNstored++;
  return true;
}

//______________________________________________________________________________
int GeantEventServer::AddEvent(GeantTaskData *td)
{
// Import next event from the generator. Thread safety has to be handled
// by the generator.
  int evt = fNload.fetch_add(1);
  if (evt >= fNevents) {
    fNload--;
    Error("AddEvent", "Event pool is full");
    return 0;
  }
  GeantEventInfo eventinfo = fRunMgr->GetPrimaryGenerator()->NextEvent();
  int ntracks = eventinfo.ntracks;
  if (!ntracks) {
    Error("AddEvent", "Event is empty");
    return 0;
  }
  fEvents[evt]->SetEvent(evt);
  fEvents[evt]->SetNprimaries(ntracks);
  fEvents[evt]->SetVertex(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);

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

  for (int itr=0; itr<ntracks; ++itr) {
    GeantTrack &track = td->GetNewTrack();
    track.fParticle = fEvents[evt]->AddPrimary(&track);
    track.SetPrimaryParticleIndex(itr); 
    track.SetPath(startpath);
    track.SetNextPath(startpath);
    track.SetEvent(evt);
    fRunMgr->GetPrimaryGenerator()->GetTrack(itr, track);
    if (!track.IsNormalized())
      track.Print("Not normalized");
    track.fBoundary = false;
    track.fStatus = kNew;
    fEvents[evt]->fNfilled++;
    if (fRunMgr->GetMCTruthMgr()) fRunMgr->GetMCTruthMgr()->AddTrack(track);
  }
  // Release path object
  VolumePath_t::ReleaseInstance(startpath);
  // Update number of stored events
  fNstored++;
  if (fNstored.load() == fNevents) {
    int nactivep = 0;
    for (int evt=0; evt<fNactiveMax; ++evt)
      nactivep += fEvents[evt]->GetNprimaries();
    // Initial basket share per propagator: nactivep/nperbasket
    fNbasketsInit = nactivep/(fRunMgr->GetConfig()->fNperBasket);
    if (!fNbasketsInit) fNbasketsInit = 1;
    if (fNbasketsInit < fRunMgr->GetNpropagators()) {
      Error("EventServer", "Too many worker threads for this configuration.");
    }
    fRunMgr->SetInitialShare(fNbasketsInit/fRunMgr->GetNpropagators());
    for (int evt=fNactiveMax; evt<fNevents; ++evt)
      nactivep += fEvents[evt]->GetNprimaries();
    fRunMgr->SetNprimaries(nactivep);
    Print("AddEvent", "Server imported %d events cumulating %d primaries", fNevents, nactivep);
    Print("EventServer", "Initial baskets to be split among propagators: %d", fNbasketsInit);
  }
  return ntracks;
}

//______________________________________________________________________________
GeantTrack *GeantEventServer::GetNextTrack()
{
// Fetch next track of the current event. Increments current event if no more
// tracks. If current event matches last activated one, resets fHasTracks flag.
// If max event fully dispatched, sets the fDone flag.

  int evt = fCurrentEvent.load();
  GeantEvent *event;
  int itr;
  while (1) {
    event = fEvents[evt];
    itr = event->fNdispatched.fetch_add(1);
    if (itr > event->GetNprimaries() - 1) {
      // Current event dispatched, try to fetch from next event
      if (evt == fLastActive.load()) {
        // No events available, check if this is last event
        fHasTracks = false;
        if (evt == fNevents-1) fDone = true;
        return nullptr;
      }
      // Attempt to change the event
      fCurrentEvent.compare_exchange_strong(evt, evt+1);
      continue;
    }
    // Track looks valid, check if it was loaded
    if (itr < event->fNfilled.load()) break;
    // Retry from the beginning
    event->fNdispatched--;
    evt = fCurrentEvent.load();
  }
  GeantTrack *track = event->GetPrimary(itr);
  track->SetEvslot(event->GetSlot());
  return track;
}

//______________________________________________________________________________
int GeantEventServer::FillBasket(GeantTrack_v &tracks, int ntracks)
{
// Fill concurrently a basket of tracks, up to the requested number of tracks.
// The client should test first the track availability using HasTracks().
  if (!fHasTracks) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack();
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
int GeantEventServer::FillBasket(Basket *basket, int ntracks)
{
// Fill concurrently a basket of tracks, up to the requested number of tracks.
// The client should test first the track availability using HasTracks().
  if (!fHasTracks) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack();
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
int GeantEventServer::FillStackBuffer(StackLikeBuffer *buffer, int ntracks)
{
// Fill concurrently up to the requested number of tracks into a stack-like buffer.
// The client should test first the track availability using HasTracks().

// *** I should template on the container to be filled, making sure that all
//     containers provide AddTrack(GeantTrack *)
  if (!fHasTracks) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack();
    if (!track) break;
    buffer->AddTrack(track);
    ndispatched++;
  }
  if (fInitialPhase) {
    int nserved = fNserved.fetch_add(1) + 1;
    if (nserved >= fNbasketsInit) fInitialPhase = false;
  }
  return ndispatched;
}

//______________________________________________________________________________
int GeantEventServer::ActivateEvents()
{
// Activate events depending on the buffer status. If the number of already
// active events is smaller than the fNactiveMax, activate as many events as
// needed to reach this, otherwise activate a single event.
  if (fEventsServed) return 0;
  int nactivated = 0;
  int nactive = fNactive.load();
  while (nactive < fNactiveMax && !fEventsServed) {
    // Try to activate next event by getting a slot
    size_t slot = 0;
    if (fFreeSlots.dequeue(slot)) {
      fNactive++;
//    if (fNactive.compare_exchange_strong(nactive, nactive+1)) {
      // We can actually activate one event
      int lastactive = fLastActive.fetch_add(1) + 1;
      fEvents[lastactive]->SetSlot(slot);
      nactivated++;
      if (lastactive == fNevents - 1) fEventsServed = true;
      fHasTracks = true;
    }
    nactive = fNactive.load();
  }
  return nactivated;
}

//______________________________________________________________________________
void GeantEventServer::CompletedEvent(int evt)
{
// Signals that event 'evt' was fully transported.
  fNactive--;
  fNcompleted++;
  if (!fFreeSlots.enqueue(fEvents[evt]->GetSlot())) Fatal("CompletedEvent", "Cannot enqueue slot");
  ActivateEvents();
}

} // GEANT_IMPL_NAMESPACE
} // Geant
