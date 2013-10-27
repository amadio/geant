#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h" 
#endif

//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
   GeantVolumeBaskets   *fManager;         // Manager for the basket
   GeantTrack_v      fTracksIn;            // Vector of input tracks
   GeantTrack_v      fTracksOut;           // Vector of output tracks
//   GeantHit_v        fHits;                // Vector of produced hits

public:
   GeantBasket();
   GeantBasket(Int_t size, GeantVolumeBaskets *mgr);
   virtual ~GeantBasket();
   
   // Add track from generator or physics process
   void              AddTrack(const GeantTrack *track);
   // Add track from a track_v array (copied)
   void              AddTrack(const GeantTrack_v &tracks, Int_t itr);
   // Add multiple tracks from a track_v array
   void              AddTracks(const GeantTrack_v &tracks, Int_t istart, Int_t iend);
   virtual void      Clear(Option_t *option="");
   Bool_t            Contains(Int_t evstart, Int_t nevents=1) const;
   Int_t             GetNinput() const {return fTracksIn.GetNtracks();}
   Int_t             GetNoutput() const {return fTracksOut.GetNtracks();}
   GeantTrack_v     &GetInputTracks() const {return (GeantTrack_v&)fTracksIn;}
   GeantTrack_v     &GetOutputTracks() const {return (GeantTrack_v&)fTracksOut;}
   Int_t             GetVolId() const {return fVolId;}
   Int_t             GetMatId() const {return fMatId;}
   void              SetVolId(Int_t id) {fVolId = id;}
   void              SetMatId(Int_t id) {fMatId = id;}
   virtual void      Print(Option_t *option="") const;
   void              PrintTrack(Int_t itr, Bool_t input=kTRUE) const;
   void              Recycle();
   
   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};   
#endif
