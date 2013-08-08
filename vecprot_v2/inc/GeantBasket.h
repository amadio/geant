#ifndef GEANT_BASKET
#define GEANT_BASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif
 
//==============================================================================
// Basket of tracks in the same volume which are transported by a single thread
//==============================================================================

class TGeoVolume;

//______________________________________________________________________________
class GeantBasket : public TObject {
protected:
   Int_t             fVolID;               // Unique TGeo volume identifier
   Int_t             fMatID;               // Material ID
   GeantTrack_v      fTracksIn;            // Vector of input tracks
   GeantTrack_v      fTracksOut;           // Vector of output tracks
//   GeantHit_v        fHits;                // Vector of produced hits

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
#endif
