#ifndef GEANT_THREADDATA
#define GEANT_THREADDATA

#ifndef ROOT_TObject
#include "TObject.h"
#endif

#ifndef GEANT_TRACK
#include "GeantTrack.h"
#endif

 
//==============================================================================
// Data organized per thread
//==============================================================================

//class TGeoHMatrix;
class TGeoVolume;
class TRandom;
class TArrayI;
class TGeoHelix;
class TGeoRotation;

//______________________________________________________________________________
class GeantThreadData : public TObject {
public:
   Int_t             fMaxPerBasket;          //! Max number of tracks per basket
   Int_t             fNprocesses;            //! Number of physics processes
   Int_t             fSizeDbl;               //! Size of array of doubles
   Int_t             fSizeBool;              //! Size of bool array
   TGeoVolume       *fVolume;                //! Current volume per thread
   TRandom          *fRndm;                  //! Random generator for thread
   Bool_t           *fBoolArray;             //![5*fMaxPerBasket] Thread array of bools
   Double_t         *fDblArray;              //![5*fMaxPerBasket] Thread array of doubles
   Double_t         *fProcStep;              //![fNprocesses*fMaxPerBasket] 
   TGeoHelix        *fFieldPropagator;       //! Propagator in magnetic field
   TGeoRotation     *fRotation;              //! Rotation used by the field propagator
   GeantTrack        fTrack;                 // Track support for this thread

public:
   GeantThreadData();
   GeantThreadData(Int_t maxperbasket, Int_t maxprocesses);
   virtual ~GeantThreadData();

   Double_t         *GetProcStep(Int_t iproc) {return fProcStep+iproc*fMaxPerBasket;}
   Double_t         *GetDblArray(Int_t size);
   Bool_t           *GetBoolArray(Int_t size);
   
private:
   GeantThreadData(const GeantThreadData &);//no imp.	
   GeantThreadData& operator=(const GeantThreadData &);//no imp.
   
   ClassDef(GeantThreadData,1)  // Stateful data organized per thread
};   
#endif
