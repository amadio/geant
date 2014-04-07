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

// class TGeoHMatrix;
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
  // TGeoHMatrix      *fMatrix;                //! Current matrix per thread
   TGeoVolume       *fVolume;                //! Current volume per thread
   TRandom          *fRndm;                  //! Random generator for thread
   Double_t         *fDblArray;              //![5*fMaxPerBasket] Support for random seed
   Double_t         *fProcStep;              //![fNprocesses*fMaxPerBasket] 
   TGeoHelix        *fFieldPropagator;       //! Propagator in magnetic field
   TGeoRotation     *fRotation;              //! Rotation used by the field propagator
   GeantTrack        fTrack;                 // Track support for this thread

public:
   GeantThreadData();
   GeantThreadData(Int_t maxperbasket, Int_t maxprocesses);
   virtual ~GeantThreadData();

   Double_t         *GetProcStep(Int_t iproc) {return fProcStep+iproc*fMaxPerBasket;}
   
   
   ClassDef(GeantThreadData,1)  // Stateful data organized per thread
};   
#endif
