#ifndef GEANT_VAPPLICATION
#define GEANT_VAPPLICATION

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantHitBlock;
class GeantTrack_v;

//______________________________________________________________________________
class GeantVApplication : public TObject {
public:
   GeantVApplication();
   virtual ~GeantVApplication() {}

   virtual void  StepManager(Int_t tid, Int_t npart, const GeantTrack_v & tracks) = 0;
   virtual void  Digitize(Int_t event) = 0;
      
   ClassDef(GeantVApplication, 1)      // User application
};
#endif
