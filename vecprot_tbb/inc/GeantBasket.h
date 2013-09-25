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
   Int_t             fNtracks;               // Number of tracks
   Int_t             fMaxTracks;             // Max number of tracks
   Int_t            *fIndex;                 //[fNtracks] Track indices in the global stack
   // have to change it to hold local tracks

public:
   GeantBasket();
   GeantBasket(Int_t maxtracks);
   virtual ~GeantBasket();

   void              AddTrack(Int_t itrack);
   void              AddTracks(const Int_t *array, Int_t ntracks);
   virtual void      Clear(Option_t *option="");
   Bool_t            Contains(Int_t event) const;
   Int_t             GetNtracks() const {return fNtracks;}
   Int_t            *GetTracks() const {return fIndex;}
   virtual void      Print(Option_t *option="") const;
   void              Resize(Int_t newSize);

   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};

//==============================================================================
// GeantTrackCollector - A class collecting all tracks propagated by a single
// thread across the boundaries of the current volume.
//==============================================================================

class GeantVolumeBasket;

//______________________________________________________________________________
class GeantTrackCollection : public TObject {
protected:
   Int_t               fNtracks;              // Number of crossing tracks
   Int_t               fSize;                 // Size of the arrays > fNtracks
   Int_t              *fTracks;               //! Indices of crossing tracks
   GeantVolumeBasket **fBaskets;              //! Volume basket to which the track goes

public:
   GeantTrackCollection();
   GeantTrackCollection(Int_t size);
   virtual ~GeantTrackCollection();

   GeantTrackCollection& operator=(const GeantTrackCollection& other);

   Int_t               AddTrack(Int_t itrack, GeantVolumeBasket *basket);
   virtual void        Clear(Option_t *option="");
   Int_t               GetNtracks() const {return fNtracks;}
   void               FlushTracks(GeantMainScheduler *main, Int_t* pushedN, Int_t* pushedP);
   virtual void        Print(Option_t *option="") const;

   ClassDef(GeantTrackCollection, 1)     // Track collection per thread
};

#endif
