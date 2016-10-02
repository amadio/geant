#include "GeantEventServer.h"

#include "globals.h"
#include "Geant/Error.h"

#include "GeantTrack.h"
#include "GeantEvent.h"
#include "GeantRunManager.h"
#include "LocalityManager.h"
#include "PrimaryGenerator.h"
#include "GeantTaskData.h"

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
  :fNevents(event_capacity), fCurrentEvent(0), fNload(0), fNstored(0), fRunMgr(runmgr)
{
// Constructor
  assert(event_capacity > 0);
  fLastActive.store(-1);
  fDoneEvents = BitSet::MakeInstance(event_capacity);
  fEvents = new GeantEvent*[event_capacity];
  for (int i=0; i<event_capacity; ++i)
    fEvents[i] = new GeantEvent();
  // Configure the locality manager if not done already
  LocalityManager *mgr = LocalityManager::Instance();
  if (!mgr->IsInitialized()) {
    mgr->SetNblocks(100);
    mgr->SetBlockSize(10000);
    mgr->SetMaxDepth(runmgr->GetConfig()->fMaxDepth);
    mgr->Init();
  }
}

//______________________________________________________________________________
GeantEventServer::~GeantEventServer()
{
// Destructor
/* Activate when numa node added to track
  LocalityManager *loc_mgr = LocalityManager::Instance();
  GeantTrack *track;
  int nprim;
*/
  for (int i = 0; i < fNload.load(); ++i) {
/* Activate when numa node added to track
    nprim = event->GetNprimaries();
    for (int j = 0; j < nprim; ++j) {
      track = fEvents[i]->GetPrimary(j);
      TrackManager &trk_mgr = loc_mgr->GetTrackManager(track->fNode);
      trk_mgr->ReleaseTrack(track);
    }
*/
    delete fEvents[i];
  }
  delete [] fEvents;
}

//______________________________________________________________________________
int GeantEventServer::AddEvent(GeantTaskData *td)
{
// Import next event from the generator. Thread safety has to be handled
// by the generator.
  int node = (td == nullptr) ? 0 : td->fNode;
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
  LocalityManager *loc_mgr = LocalityManager::Instance();
  TrackManager &trk_mgr = loc_mgr->GetTrackManager(node);

  Volume_t *vol = 0;
  // Initialize the start path
  VolumePath_t *startpath = VolumePath_t::MakeInstance(fRunMgr->GetConfig()->fMaxDepth);
#ifdef USE_VECGEOM_NAVIGATOR
  vecgeom::SimpleNavigator nav;
  startpath->Clear();
  nav.LocatePoint(GeoManager::Instance().GetWorld(),
                  Vector3D<Precision>(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert), *startpath, true);
  vol = const_cast<Volume_t *>(startpath->Top()->GetLogicalVolume());
#else
  TGeoNavigator *nav = gGeoManager->GetCurrentNavigator();
  if (!nav)
    nav = gGeoManager->AddNavigator();
  TGeoNode *node = nav->FindNode(eventinfo.xvert, eventinfo.yvert, eventinfo.zvert);
  vol = node->GetVolume();
  startpath->InitFromNavigator(nav);
#endif
  if (td) td->fVolume = vol;
  
  for (int itr=0; itr<ntracks; ++itr) {
    GeantTrack &track = trk_mgr.GetTrack();
    fEvents[evt]->AddPrimary(&track);
    track.SetPath(startpath);
    track.SetNextPath(startpath);
    track.SetEvent(evt);
    fRunMgr->GetPrimaryGenerator()->GetTrack(itr, track);
    if (!track.IsNormalized())
      track.Print("Not normalized");
    track.fBoundary = false;
    track.fStatus = kAlive;
    fEvents[evt]->fNfilled++;
  }
  // Release path object
  VolumePath_t::ReleaseInstance(startpath);
  // Update number of stored events
  fNstored++;
  return ntracks;
}

//______________________________________________________________________________
GeantTrack *GeantEventServer::GetNextTrack()
{
// Fetch next track
  int evt = fCurrentEvent.load();
  GeantEvent *event;
  int itr;
  while (1) {
    event = fEvents[evt];
    itr = event->fNdispatched.fetch_add(1);
    if (itr > event->GetNprimaries() - 1) {
      if (evt == fNevents-1) {
        // If last event just exit
        fDone = true;
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
  return event->GetPrimary(itr);
}

//______________________________________________________________________________
int GeantEventServer::FillBasket(GeantTrack_v &tracks, int ntracks)
{
// Fill concurrently a basket of tracks
  if (fDone) return 0;
  int ndispatched = 0;
  for (int i=0; i<ntracks; ++i) {
    GeantTrack *track = GetNextTrack();
    if (!track) break;
    tracks.AddTrack(*track);
    ndispatched++;
  }
  return ndispatched;
}

//______________________________________________________________________________
bool GeantEventServer::ActivateEvents()
{
// Activate events depending on the buffer status. If the number of already
// active events is smaller than the fNactiveMax, activate as many events as
// needed to reach this, otherwise activate a single event.
  
}

} // GEANT_IMPL_NAMESPACE
} // Geant
