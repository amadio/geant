#ifndef GXTrackHandler_H
#define GXTrackHandler_H 1

#include "GXTrack.h"

class GXTrackHandler 
{
public:
 
  GXTrackHandler();
  GXTrackHandler(size_t nTracks, double photonFrac = 0.0);
  ~GXTrackHandler();

  void   SetNumberOfTracks(size_t nTracks, double photonFrac = 0.0);
  size_t GetNumberOfTracks() {return theNumberOfTracks; };
  size_t GetNumberOfElectrons() {return theNumberOfElectrons; };
  size_t GetNumberOfPhotons() {return theNumberOfPhotons; };

  void Allocate(size_t nTracks, double photonFrac = 0.0);
  void Deallocate();
  void Reallocate(size_t nTracks, double photonFrac = 0.0);

  GXTrack* GetTracks() {return track_h; };
  GXTrack* GetElectrons() {return track_e; };
  GXTrack* GetPhotons() {return track_g; };

  void FillOneTrack(GXTrack* aTrack);
  void GenerateRandomTracks(size_t nTracks, double photonFrac = 0.0,
			    double minP=20.0, double maxP = 1000.);
  void CopyTrack(GXTrack *This, GXTrack *track);

private:

  size_t theNumberOfTracks;
  size_t theNumberOfElectrons;
  size_t theNumberOfPhotons;

  double thePhotonFraction;

  GXTrack *track_h;
  GXTrack *track_e;
  GXTrack *track_g;

};
#endif
