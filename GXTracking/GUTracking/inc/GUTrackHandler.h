#ifndef GUTrackHandler_H
#define GUTrackHandler_H 1

#include "GUTrack.h"
#include "stddef.h"

#include "backend/Backend.h"


namespace vecphys {

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

  GUTrack*   GetAoSTracks() {return fTrack_aos; };
  GUTrack_v& GetSoATracks() {return fTrack_soa; };

  void FillOneTrack(GUTrack* aTrack);
  void GenerateRandomTracks(size_t nTracks, 
			    double minP=20.0, double maxP = 1000.);
  void CopyAoSTracks(GUTrack *fromAoS, GUTrack *toAoS);
  void CopySoATracks(GUTrack_v& fromSoA, GUTrack_v& toSoA);

private:

  size_t fNumberOfTracks;
  GUTrack*  fTrack_aos;
  GUTrack_v fTrack_soa;
  char*     fBuffer;
};

} // end namespace vecphys

#endif
