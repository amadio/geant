#ifndef GEANT_VAPPLICATION
#define GEANT_VAPPLICATION

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantHitBlock;

//______________________________________________________________________________
class GeantVApplication : public TObject {
public:
   GeantVApplication();
   virtual ~GeantVApplication() {}

   virtual void  StepManager(Int_t tid, Int_t iproc, Int_t npart, Int_t *particles) = 0;
   virtual void  Digitize(Int_t event) = 0;
      
   ClassDef(GeantVApplication, 1)      // User application
};
#endif
