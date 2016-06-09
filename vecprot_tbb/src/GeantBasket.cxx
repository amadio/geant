#include "GeantBasket.h"

#include "GeantVolumeBasket.h"
#include "GeantMainScheduler.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"

#include "TThread.h"

//______________________________________________________________________________
GeantBasket::GeantBasket()
  : TObject(), fNtracks(0), fMaxTracks(0), fIndex(0) {
  // ctor.
}

//______________________________________________________________________________
GeantBasket::GeantBasket(int maxtracks) : TObject(), fNtracks(0), fMaxTracks(maxtracks), fIndex(new int[maxtracks]) {
  // ctor.
}

//______________________________________________________________________________
GeantBasket::~GeantBasket() {
  // dtor.
  delete[] fIndex;
}

//______________________________________________________________________________
void GeantBasket::AddTrack(int itrack) {
  // Add a new track to this basket;
  if (fNtracks >= fMaxTracks - 1)
    Resize(2 * fMaxTracks);
  fIndex[fNtracks++] = itrack;
}

//______________________________________________________________________________
void GeantBasket::AddTracks(const int *array, int ntracks) {
  // Add array of tracks to the basket.
  if (fNtracks + ntracks > fMaxTracks - 1)
    Resize(TMath::Max(fNtracks + ntracks, 2 * fMaxTracks));
  memcpy(&fIndex[fNtracks], array, ntracks * sizeof(int));
  fNtracks += ntracks;
}

//______________________________________________________________________________
void GeantBasket::Clear(const char *) {
  // Clear basket;
  fNtracks = 0;
}

//______________________________________________________________________________
bool GeantBasket::Contains(int event) const {
  // Checks if any of the array of tracks belongs to the given event.

  GeantPropagator *gPropagator = GeantPropagator::Instance();

  for (int itr = 0; itr < fNtracks; itr++) {
    if (gPropagator->fTracks[fIndex[itr]]->event == event)
      return kTRUE;
  }
  return kFALSE;
}

//______________________________________________________________________________
void GeantBasket::Print(const char *) const {
  // Print basket content.
  TThread::Lock();
  TString s = Form("basket %p with %d tracks:", this, fNtracks);
  for (int i = 0; i < fNtracks; i++) {
    //      gPropagator->fTracks[fIndex[i]]->Print();
    s += Form(" %d", fIndex[i]);
  }
  Printf("%s", s.Data());
  TThread::UnLock();
}

//______________________________________________________________________________
void GeantBasket::Resize(int newSize) {
  // Resize the array of track indices. Not thread safe - should be called by a
  // single thread at a time;
  int *newindex = new int[newSize];
  memcpy(newindex, fIndex, fNtracks * sizeof(int));
  delete[] fIndex;
  fIndex = newindex;
  fMaxTracks = newSize;
}

//______________________________________________________________________________
GeantTrackCollection::GeantTrackCollection()
  : TObject(), fNtracks(0), fSize(0), fTracks(0), fBaskets(0) {
  // Default ctor.
}

//______________________________________________________________________________
GeantTrackCollection::GeantTrackCollection(int size) : TObject(), fNtracks(0), fSize(size), fTracks(0), fBaskets(0) {
  // Ctor.
  fTracks = new int[size];
  fBaskets = new GeantVolumeBasket *[size];
}

//______________________________________________________________________________
GeantTrackCollection::~GeantTrackCollection() {
  // Dtor.
  delete[] fTracks;
  delete[] fBaskets;
}

GeantTrackCollection &GeantTrackCollection::operator=(const GeantTrackCollection &other) {
  if (&other != this) {
    fNtracks = other.fNtracks;
    fSize = other.fSize;

    if (fTracks)
      delete[] fTracks;
    fTracks = new int[fSize];
    memcpy(fTracks, other.fTracks, fSize * sizeof(int));

    if (fBaskets)
      delete[] fBaskets;
    fBaskets = new GeantVolumeBasket *[fSize];
    memcpy(fBaskets, other.fBaskets, fSize * sizeof(GeantVolumeBasket *));
  }
  return *this;
}

//______________________________________________________________________________
void GeantTrackCollection::Clear(const char *) {
  // Clear basket;
  fNtracks = 0;
}

//______________________________________________________________________________
int GeantTrackCollection::AddTrack(int itrack, GeantVolumeBasket *basket) {
  GeantPropagator *propagator = GeantPropagator::Instance();

  // Add a new track entering the basket.
  if (!propagator->fTracks[itrack]->IsAlive())
    return fNtracks;
  if (fNtracks == fSize - 1) {
    int *tracks = new int[2 * fSize];
    GeantVolumeBasket **baskets = new GeantVolumeBasket *[2 * fSize];
    memcpy(tracks, fTracks, fSize * sizeof(int));
    memcpy(baskets, fBaskets, fSize * sizeof(GeantVolumeBasket *));
    delete[] fTracks;
    fTracks = tracks;
    delete[] fBaskets;
    fBaskets = baskets;
    fSize *= 2;
  }
  fTracks[fNtracks] = itrack;
  fBaskets[fNtracks] = basket;
  return fNtracks++;
}

//______________________________________________________________________________
void GeantTrackCollection::FlushTracks(GeantMainScheduler *main, int *pushedN, int *pushedP) {
  // Flush all tracks to the main scheduler. Returns number of injected baskets.
  bool prior;
  int injected;
  for (int itr = 0; itr < fNtracks; itr++) {
    injected = main->AddTrack(fTracks[itr], fBaskets[itr]->GetNumber(), &prior);
    if (injected) {
      if (!prior)
        (*pushedN)++;
      else
        (*pushedP)++;
    }
  }
  fNtracks = 0;
  return;
}

//______________________________________________________________________________
void GeantTrackCollection::Print(const char *) const {
  // Print info
  TThread::Lock();
  TString s = Form("collection %p with %d tracks:", this, fNtracks);
  for (int i = 0; i < fNtracks; i++) {
    s += Form(" %d", fTracks[i]);
  }
  Printf("%s", s.Data());
  TThread::UnLock();
}
