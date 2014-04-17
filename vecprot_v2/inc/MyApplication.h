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

class GeantTrack_v;

//______________________________________________________________________________
class MyApplication : public GeantVApplication {
   static const Int_t kNlayers = 10;
private:
   Bool_t               fInitialized; // initialized flag
   Int_t                fIdGap;       // ID for the gap volume
   Int_t                fIdAbs;       // ID for the absorber volume
   Float_t              fEdepGap[kNlayers]; // Energy deposition per layer
   Float_t              fLengthGap[kNlayers]; // step length in every layer
   Float_t              fEdepAbs[kNlayers]; // Energy deposition per layer
   Float_t              fLengthAbs[kNlayers]; // step length in every layer
   GeantFactory<MyHit> *fFactory; // Hits factory      

public:
   MyApplication();
   virtual ~MyApplication() {}

   virtual Bool_t Initialize();
   virtual void   StepManager(Int_t tid, Int_t npart, const GeantTrack_v & tracks);
   virtual void   Digitize(Int_t event);
      
   ClassDef(MyApplication, 1)      // User application
};
#endif
