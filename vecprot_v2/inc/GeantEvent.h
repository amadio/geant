#ifndef GEANT_EVENT
#define GEANT_EVENT

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

class GeantEvent : public TObject {
private:
   Int_t  fEvent;        // event number
   Int_t  fSlot;         // fixed slot number
   Int_t  fNtracks;      // number of tracks
   Int_t  fNdone;        // number of done tracks
   TMutex fMutex;        // mutex for this event

public:
   GeantEvent() : fEvent(0), fSlot(0), fNtracks(0), fNdone(0), fMutex() {}
   GeantEvent(Int_t ntr) : fEvent(0), fSlot(0), fNtracks(ntr), fNdone(0), fMutex() {}
   ~GeantEvent() {}
   
   Int_t              AddTrack();
   Int_t              GetEvent() const {return fEvent;}
   Int_t              GetSlot() const {return fSlot;}
   Int_t              GetNtracks() const {return fNtracks;}
   void               SetEvent(Int_t event) {fEvent = event;}
   void               SetSlot(Int_t islot) {fSlot = islot;}
   void               Reset() {fNtracks = fNdone = 0;}
   void               StopTrack();
   Bool_t             Transported() const {return ((fNtracks>0) && (fNtracks==fNdone));}

   ClassDef(GeantEvent, 1)       // The G5 event
};
#endif
