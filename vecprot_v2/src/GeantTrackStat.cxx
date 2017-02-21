#include "GeantTrackStat.h"
#include "GeantTrackVec.h"
#include "Geant/Error.h"

//______________________________________________________________________________
GeantTrackStat::GeantTrackStat(int nslots)
  : fNslots(nslots), fNtracks(0), fNsteps(0), fMutex() {
  // Ctor
  InitArrays(nslots);
}

//______________________________________________________________________________
GeantTrackStat::~GeantTrackStat() {
  // Dtor.
  delete[] fNtracks;
  delete[] fNsteps;
}

//______________________________________________________________________________
void GeantTrackStat::AddTrack(const GeantTrack &track) {
  // Add a track
  fNtracks[track.fEvslot]++;
  fNsteps[track.fEvslot] += track.fNsteps;
}

//______________________________________________________________________________
void GeantTrackStat::AddTrack(const GeantTrack_v &trackv, int itr) {
  // Add a track from an array
  fNtracks[trackv.fEvslotV[itr]]++;
  fNsteps[trackv.fEvslotV[itr]] += trackv.fNstepsV[itr];
}

//______________________________________________________________________________
void GeantTrackStat::AddTracks(const GeantTrack_v &trackv) {
  // Remove statistics for tracks
  fMutex.lock();
  int ntracks = trackv.GetNtracks();
  for (int i = 0; i < ntracks; i++) {
    fNtracks[trackv.fEvslotV[i]]++;
    fNsteps[trackv.fEvslotV[i]] += trackv.fNstepsV[i];
  }
  fMutex.unlock();
}

//______________________________________________________________________________
void GeantTrackStat::RemoveTracks(const GeantTrack_v &trackv) {
  // Remove statistics for tracks
  fMutex.lock();
  int ntracks = trackv.GetNtracks();
  for (int i = 0; i < ntracks; i++) {
    // do *NOT* remove new tracks since they were not added yet anywhere
    if (trackv.fStatusV[i] == Geant::kNew)
      continue;
    fNtracks[trackv.fEvslotV[i]]--;
    fNsteps[trackv.fEvslotV[i]] -= trackv.fNstepsV[i];
  }
  fMutex.unlock();
}

//______________________________________________________________________________
GeantTrackStat &GeantTrackStat::operator+=(const GeantTrackStat &other) {
  // Compound addition.
  if (fNslots != other.fNslots) {
    Geant::Error("GeantTrackStat::operator+=", "Different number of slots");
    return *this;
  }
  for (int i = 0; i < fNslots; i++) {
    fNtracks[i] += other.fNtracks[i];
    fNsteps[i] += other.fNsteps[i];
  }
  return *this;
}

//______________________________________________________________________________
GeantTrackStat &GeantTrackStat::operator-=(const GeantTrackStat &other) {
  // Compound addition.
  if (fNslots != other.fNslots) {
    Geant::Error("GeantTrackStat::operator+=", "Different number of slots");
    return *this;
  }
  for (int i = 0; i < fNslots; i++) {
    fNtracks[i] -= other.fNtracks[i];
    fNsteps[i] -= other.fNsteps[i];
  }
  return *this;
}

//______________________________________________________________________________
void GeantTrackStat::InitArrays(int nslots) {
  // Initialize arrays
  fNslots = nslots;
  delete[] fNtracks;
  delete[] fNsteps;
  fNtracks = new int[nslots];
  memset(fNtracks, 0, nslots * sizeof(int));
  fNsteps = new int[nslots];
  memset(fNsteps, 0, nslots * sizeof(int));
}

//______________________________________________________________________________
void GeantTrackStat::Print(const char *) const {
  // Print statistics
  printf("slot: ");
  for (int i = 0; i < fNslots; i++) {
    printf("%5d ", i);
  }
  printf("\n");
  printf("ntr:  ");
  for (int i = 0; i < fNslots; i++) {
    printf("%5d ", fNtracks[i]);
  }
  printf("\n");
}

//______________________________________________________________________________
void GeantTrackStat::Reset() {
  // Reset statistics
  if (fNslots) {
    memset(fNtracks, 0, fNslots * sizeof(int));
    memset(fNsteps, 0, fNslots * sizeof(int));
  }
}
