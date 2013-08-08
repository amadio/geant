#ifndef GEANT_HITBLOCK
#define GEANT_HITBLOCK

#ifndef ROOT_TObject
#include "TObject.h"
#endif

//______________________________________________________________________________
// GeantHitBlock - a fixed-size block of user hits stored contiguously in memory
//                 The class is used to efficiently manage hit user data in a
//                 multithreaded environment.
//______________________________________________________________________________

class TClass;

class GeantHitBlock : public TObject {
public:
   Int_t                 fSize;      // fixed size
   Int_t                 fNext;      // number of hits
   TClass               *fClass;     // class type for the hits (MyHit::Class())
   Int_t                 fClSize;    // size of hit class stored
   char                 *fSupport;   // support array
   void                 *fArray;     // hits array
   
   GeantHitBlock(TClass *cl, Int_t size, void *storage=0);
   virtual ~GeantHitBlock();
   
   Bool_t        Full() const {return (fNext==fSize);}
   TObject      *NextHit();
   void          ClearHits();
   void          CopyHit(Int_t ihit, GeantHitBlock &other);
   
   ClassDef(GeantHitBlock, 0)  // A block of hits
};
#endif
