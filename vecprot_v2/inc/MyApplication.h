#ifndef GEANT_MYAPPLICATION
#define GEANT_MYAPPLICATION

#ifndef GEANT_VAPPLICATION
#include "GeantVApplication.h"
#endif

#ifndef GEANT_FACTORY
#include "GeantFactory.h"
#endif

#ifndef GEANT_MYHIT
#include "MyHit.h"
#endif

//______________________________________________________________________________
class MyApplication : public GeantVApplication {
private:
   GeantFactory<MyHit> *fFactory; // Hits factory      

public:
   MyApplication();
   virtual ~MyApplication() {}

   virtual void  StepManager(Int_t tid, Int_t iproc, Int_t npart, Int_t *particles);
   virtual void  Digitize(Int_t event);
      
   ClassDef(MyApplication, 1)      // User application
};
#endif
