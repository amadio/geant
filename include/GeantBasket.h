#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include "TMutex.h"
 
//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

class TGeoVolume;

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
   Int_t             fNtracks;               // Number of tracks
   Int_t             fMaxTracks;             // Max number of tracks
   Int_t            *fIndex;                 //[fNtracks] Track indices in the global stack

public:
   GeantBasket();
   GeantBasket(Int_t maxtracks);
   virtual ~GeantBasket();
   
   void              AddTrack(Int_t itrack);
   void              AddTracks(const Int_t *array, Int_t ntracks);
   virtual void      Clear(Option_t *option="");
   Int_t             GetNtracks() const {return fNtracks;}
   Int_t            *GetTracks() const {return fIndex;}
   virtual void      Print(Option_t *option="") const;
   void              Resize(Int_t newSize);
   
   ClassDef(GeantBasket,1)  // A basket containing tracks in the same geomety volume
};   

//==============================================================================
// Class holding a list of track baskets that can be used in a preferential
// round robin scheduling scheme.
//==============================================================================

//______________________________________________________________________________
class GeantBasketScheduler : public TObject {
protected:
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNbaskets;              // Number of baskets
   Int_t             fMaxBaskets;            // Maximum number of baskets
   Int_t             fVIB;                   // Very important basket getting the priority
   Int_t             fCrt;                   // Current basket to be filled
   Int_t             fNbusy;                 // Number of busy baskets
   GeantBasket     **fBaskets;               //[fNbaskets] Array of baskets
   GeantBasket      *fCollector;             //! garbage collector basket
   TMutex            fBasketLock;            //! Mutex lock for the basket

public:
   GeantBasketScheduler();
   GeantBasketScheduler(TGeoVolume *vol);
   virtual ~GeantBasketScheduler();

   void              AddTrack(Int_t itrack);
   void              AddNewBasket();
   GeantBasket      *GarbageCollect(GeantBasket *basket, Bool_t feed);
   Int_t             GetNtracks() const;
   Int_t             GetNtotal() const;
   TGeoVolume       *GetVolume() const {return fVolume;}
   void              InjectBasket(Int_t islot);

   ClassDef(GeantBasketScheduler,1)  // A basket scheduler
};
      
#endif
