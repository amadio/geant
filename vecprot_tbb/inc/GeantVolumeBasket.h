#ifndef GEANT_VOLUMEBASKET
#define GEANT_VOLUMEBASKET

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef ROOT_TGeoVolume
#include "TGeoVolume.h"
#endif

#ifndef GEANT_BASKET
#include "GeantBasket.h"
#endif

#include "GeantTrack.h"

class TGeoBranchArray;

//______________________________________________________________________________
class GeantVolumeBasket : public TObject {
protected:
   TGeoVolume       *fVolume;                // Volume for which applies
   Int_t             fNumber;                // Number assigned

public:
   GeantVolumeBasket(TGeoVolume *vol, Int_t number);
   virtual ~GeantVolumeBasket();

   const char       *GetName() const              {return (fVolume)?fVolume->GetName():ClassName();}
   Int_t             GetNumber() const            {return fNumber;}
   TGeoVolume       *GetVolume() const            {return fVolume;}
   virtual void      Print(Option_t *option="") const;

   void              ComputeTransportLength(Int_t ntracks, Int_t *trackin);
   void              PropagateTracks(Int_t ntracks, Int_t *trackin, Int_t &nout, Int_t *trackout, Int_t &ntodo, Int_t *tracktodo, Int_t &ncross, Int_t *trackcross);
   static void       ResetStep(Int_t ntracks, Int_t *array);

   ClassDef(GeantVolumeBasket,1)  // A path in geometry represented by the array of indices
};

#endif // GEANT_VOLUMEBASKET

