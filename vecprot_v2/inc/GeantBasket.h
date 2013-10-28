#ifndef GEANT_BASKET
#define GEANT_BASKET

#include "TObject.h"
#include "TGeoExtension.h"
#include "TGeoVolume.h"
#include "GeantTrack.h" 
#include "sync_objects.h"

//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

class TGeoVolume;
class GeantBasketMgr;

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
   GeantBasketMgr   *fManager;         // Manager for the basket
   GeantTrack_v      fTracksIn;            // Vector of input tracks
   GeantTrack_v      fTracksOut;           // Vector of output tracks
//   GeantHit_v        fHits;                // Vector of produced hits

public:
   GeantBasket();
   GeantBasket(Int_t size, GeantBasketMgr *mgr);
   virtual ~GeantBasket();
   
   // Add track from generator or physics process
   void              AddTrack(const GeantTrack &track);
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
   GeantBasketMgr   *GetBasketMgr() const {return fManager;}
   TGeoVolume       *GetVolume() const;
   virtual void      Print(Option_t *option="") const;
   void              PrintTrack(Int_t itr, Bool_t input=kTRUE) const;
   void              Recycle();
   
   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};   

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

//______________________________________________________________________________
class GeantBasketMgr : public TGeoExtension {
protected:
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNumber;                // Number assigned
   Int_t             fThreshold;             // Adjustable transportability threshold
   Int_t             fNbaskets;              // Number of baskets for this volume
   Int_t             fNused;                 // Number of baskets in use
   GeantBasket      *fCBasket;               // Current basket being filled
   GeantBasket      *fPBasket;               // Current priority basket being filled
   dcqueue<GeantBasket> fBaskets;            // queue of available baskets
   dcqueue<GeantBasket> *fFeeder;            // feeder queue to which baskets get injected
   TMutex            fMutex;                 // Mutex for this basket manager

   GeantBasket      *GetNextBasket();

public:
   GeantBasketMgr() : fVolume(0), fNumber(0), fThreshold(0), fNbaskets(0), 
                         fNused(0), fCBasket(0), fPBasket(0), fBaskets(), fFeeder(0), fMutex() {}
   GeantBasketMgr(TGeoVolume *vol, Int_t number);
   virtual ~GeantBasketMgr();
   
   virtual TGeoExtension *Grab() {return 0;}
   virtual void           Release() const {}
   Int_t             AddTrack(const GeantTrack &track, Bool_t priority=kFALSE);
   Int_t             AddTrack(const GeantTrack_v &trackv, Int_t itr, Bool_t priority=kFALSE);
   Int_t             CollectPrioritizedTracks(Int_t evmin, Int_t evmax);
   Int_t             FlushPriorityBasket();
   Int_t             GarbageCollect();
   Int_t             GetNbaskets() const          {return fNbaskets;}
   Int_t             GetNused() const             {return fNused;}
   Int_t             GetThreshold() const         {return fThreshold;}
   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNumber() const            {return fNumber;}
   TGeoVolume       *GetVolume() const            {return fVolume;}
   virtual void      Print(Option_t *option="") const;
   void              RecycleBasket(GeantBasket *b);
   void              SetFeederQueue(dcqueue<GeantBasket> *queue) {fFeeder = queue;}
   void              SetThreshold(Int_t thr)      {fThreshold = thr;}
   
   ClassDef(GeantBasketMgr,1)  // A path in geometry represented by the array of indices
};   
#endif
