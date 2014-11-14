#ifndef GUTrackHandler_H
#define GUTrackHandler_H 1

#include "GUTrack.h"
#include "stddef.h"

class GUTrackHandler 
{
public:
 
  GUTrackHandler();
  GUTrackHandler(size_t nTracks);
  ~GUTrackHandler();

  void   SetNumberOfTracks(size_t nTracks);
  size_t GetNumberOfTracks() {return fNumberOfTracks; };

  void Allocate(size_t nTracks);
  void Deallocate();
  void Reallocate(size_t nTracks);

  GUTrack* GetTracks() {return fTrack_aos; };

  void FillOneTrack(GUTrack* aTrack);
  void GenerateRandomTracks(size_t nTracks, 
			    double minP=20.0, double maxP = 1000.);
  void CopyTrack(GUTrack *This, GUTrack *track);

private:

  size_t fNumberOfTracks;
  GUTrack *fTrack_aos;
  GUTrack_v *fTrack_soa;
};

#endif
