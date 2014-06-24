#ifndef GEANT_EVENT
#define GEANT_EVENT

#if __cplusplus >= 201103L
#include <atomic>
#endif

#ifndef ROOT_TMutex
#include "TMutex.h"
#endif

class GeantEvent : public TObject {
private:
#if __cplusplus >= 201103L
   std::atomic_int  fEvent;        // event number
   std::atomic_int  fSlot;         // fixed slot number
   std::atomic_int  fNtracks;      // number of tracks
   std::atomic_int  fNdone;        // number of done tracks
#else
// dummy declarations to cheat CINT
   int  fEvent;        // event number
   int  fSlot;         // fixed slot number
   int  fNtracks;      // number of tracks
   int  fNdone;        // number of done tracks
#endif      
   TMutex fMutex;        // mutex for this event

public:
   GeantEvent() : fEvent(0), fSlot(0), fNtracks(0), fNdone(0), fMutex() {}
   GeantEvent(Int_t ntr) : fEvent(0), fSlot(0), fNtracks(ntr), fNdone(0), fMutex() {}
   ~GeantEvent() {}
   
   Int_t              AddTrack();
#if __cplusplus >= 201103L
   Int_t              GetEvent() const {return fEvent.load();}
   Int_t              GetSlot() const {return fSlot.load();}
   Int_t              GetNtracks() const {return fNtracks.load();}
   void               SetEvent(Int_t event) {fEvent.store(event);}
   void               SetSlot(Int_t islot) {fSlot.store(islot);}
   void               Reset() {fNtracks.store(0); fNdone.store(0);}
   Bool_t             Transported() const {return ((fNtracks.load()>0) && (fNtracks==fNdone));}
   void               StopTrack() {fNdone++;}
#else
   Int_t              GetEvent() const {return fEvent;}
   Int_t              GetSlot() const {return fSlot;}
   Int_t              GetNtracks() const {return fNtracks;}
   void               SetEvent(Int_t event) {fEvent = event;}
   void               SetSlot(Int_t islot) {fSlot = islot;}
   void               Reset() {fNtracks = fNdone = 0;}
   Bool_t             Transported() const {return ((fNtracks>0) && (fNtracks==fNdone));}
   void               StopTrack() {fMutex.Lock(); fNdone++; fMutex.UnLock();}
#endif
   void               Print(Option_t *option="") const;

   ClassDef(GeantEvent, 1)       // The G5 event
};
#endif
