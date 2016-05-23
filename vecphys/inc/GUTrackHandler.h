#ifndef GUTrackHandler_H
#define GUTrackHandler_H 1

#include "GUTrack.h"
#include "stddef.h"

#include "base/VPGlobal.h"

namespace vecphys {

class GUTrackHandler {
public:
  GUTrackHandler();
  GUTrackHandler(size_t nTracks);
  ~GUTrackHandler();

  void SetNumberOfTracks(size_t nTracks); // Dangerous - must be *private*
  size_t GetNumberOfTracks() { return fNumberOfTracks; };

  void Allocate(size_t nTracks);
  void Deallocate();
  void Reallocate(size_t nTracks);

  GUTrack *GetAoSTracks() { return fTrack_aos; };
  GUTrack_v &GetSoATracks() { return fTrack_soa; };

  void FillOneTrack(GUTrack *aTrack);
  void GenerateRandomTracks(size_t nTracks, double minP = 20.0, double maxP = 1000.);

  // utility functions - can be elsewhere
  void SortAoSTracksByEnergy(GUTrack *AoS, size_t nTracks);
  void SortSoATracksByEnergy(GUTrack_v &SoA, size_t nTracks);

  void CopyAoSTracks(GUTrack *fromAoS, GUTrack *toAoS, size_t Tracks);
  void CopySoATracks(GUTrack_v &fromSoA, GUTrack_v &toSoA, size_t nTracks);
  void CopyAoSTracksToSoA(GUTrack *fromAoS, GUTrack_v &toSoA, size_t nTracks);
  void CopySoATracksToAoS(GUTrack_v &fromSoA, GUTrack *toAoS, size_t nTracks);

private:
  size_t fNumberOfTracks;
  GUTrack *fTrack_aos;
  GUTrack_v fTrack_soa;
  char *fBuffer;
};

} // end namespace vecphys

#endif
