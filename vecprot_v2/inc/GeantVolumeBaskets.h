#ifndef GEANT_VOLUMEBASKETS
#define GEANT_VOLUMEBASKETS

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include "GeantTrack.h"
 
#ifndef ROOT_TGeoVolume
#include "TGeoVolume.h"
#endif

#ifndef GEANT_BASKET
#include "GeantBasket.h"
#endif

#include "sync_objects.h"

//______________________________________________________________________________
// Basket manager for a given volume. Holds a list of free baskets stored in a
// concurrent queue
//______________________________________________________________________________

//______________________________________________________________________________
class GeantVolumeBaskets : public TObject {
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
   GeantVolumeBaskets() {fVolume(0), fNumber(0), fThreshold(0), fNbaskets(0), 
                         fNused(0), fCBasket(0), fPBasket(0), fBaskets(), fFeeder(0), fMutex() {}
   GeantVolumeBaskets(TGeoVolume *vol, Int_t number);
   virtual ~GeantVolumeBasket();
   
   Int_t             AddTrack(const GeantTrack &track, Bool_t priority=kFALSE);
   Int_t             GetNbaskets() const          {return fNbaskets;}
   Int_t             GetNused() const             {return fNused;}
   Int_t             GetThreshold() const         {return fThreshold;}
   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNumber() const            {return fNumber;}
   TGeoVolume       *GetVolume() const            {return fVolume;}
   virtual void      Print(Option_t *option="") const;
   void              RecycleBasket(GeantBasket *b);
   void              SetFeederQueue(dcqueue<GeantBasket> *queue) {fFeeder = queue;}
   void              SetThreshold(Int_t thr)      {fThreshold = threshold;}
   
   ClassDef(GeantVolumeBasket,1)  // A path in geometry represented by the array of indices
};   
#endif
