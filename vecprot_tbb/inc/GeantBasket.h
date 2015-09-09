#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

/*class TGeoVolume;*/
class GeantMainScheduler;

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
  int fNtracks;   // Number of tracks
  int fMaxTracks; // Max number of tracks
  int *fIndex;    //[fNtracks] Track indices in the global stack
                  // have to change it to hold local tracks

public:
  GeantBasket();
  GeantBasket(int maxtracks);
  virtual ~GeantBasket();

  void AddTrack(int itrack);
  void AddTracks(const int *array, int ntracks);
  virtual void Clear(const char *option = "");
  bool Contains(int event) const;
  int GetNtracks() const { return fNtracks; }
  int *GetTracks() const { return fIndex; }
  virtual void Print(const char *option = "") const;
  void Resize(int newSize);

  ClassDef(GeantBasket, 1) // A basket containing tracks in the same geomety volume
};

//==============================================================================
// GeantTrackCollector - A class collecting all tracks propagated by a single
// thread across the boundaries of the current volume.
//==============================================================================

class GeantVolumeBasket;

//______________________________________________________________________________
class GeantTrackCollection : public TObject {
protected:
  int fNtracks;                 // Number of crossing tracks
  int fSize;                    // Size of the arrays > fNtracks
  int *fTracks;                 //! Indices of crossing tracks
  GeantVolumeBasket **fBaskets; //! Volume basket to which the track goes

public:
  GeantTrackCollection();
  GeantTrackCollection(int size);
  virtual ~GeantTrackCollection();

  GeantTrackCollection &operator=(const GeantTrackCollection &other);

  int AddTrack(int itrack, GeantVolumeBasket *basket);
  virtual void Clear(const char *option = "");
  int GetNtracks() const { return fNtracks; }
  void FlushTracks(GeantMainScheduler *main, int *pushedN, int *pushedP);
  virtual void Print(const char *option = "") const;

  ClassDef(GeantTrackCollection, 1) // Track collection per thread
};

#endif
