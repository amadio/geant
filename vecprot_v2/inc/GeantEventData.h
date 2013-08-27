#ifndef GEANT_THREADDATA
#define GEANT_THREADDATA

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#include <vector>
 
//==============================================================================
// Data organized per thread
//==============================================================================

class GeantEvent;
class GeantHitBlock;

//______________________________________________________________________________
class GeantEventData : public TObject {
public:
   GeantEvent              *fEvent;            // Event pointer
   vector<GeantHitBlock *>  fHitBlocks;        //! Hit blocks
// GeantDigitBlock   **fDigitBlocks;           //! Digit blocks   
public:
   GeantEventData();
   virtual ~GeantEventData();
   
   void     AddHitBlock(GeantHitBlock *block);
   void     ClearHits();
   
   ClassDef(GeantEventData,1)  // Data per event
};   
#endif
