#ifndef GEANT_THREADDATA
#define GEANT_THREADDATA

#include "Rtypes.h"

#include <vector>

//==============================================================================
// Data organized per thread
//==============================================================================

class TGeoHMatrix;
class TGeoVolume;
class TRandom;
class TArrayI;
class TGeoHelix;
class TGeoRotation;
class GeantTrackCollection;

//______________________________________________________________________________
class GeantThreadData
{
public:
   Int_t             fMaxPerBasket;          //! Max number of tracks per basket
   Int_t             fNprocesses;            //! Number of physics processes
   TGeoHMatrix      *fMatrix;                //! Current matrix per thread
   TGeoVolume       *fVolume;                //! Current volume per thread
   TRandom          *fRndm;                  //! Random generator for thread
   Double_t         *fDblArray;              //![5*fMaxPerBasket] Support for random seed
   Double_t         *fProcStep;              //![fNprocesses*fMaxPerBasket]
   TArrayI          *fPartInd;               //![fMaxPerBasket] Indices array for transported particles
   TArrayI          *fPartNext;              //![fMaxPerBasket] Indices array for remaining particles
   TArrayI          *fPartTodo;              //![fMaxPerBasket] Indices array for todo particles
   TArrayI          *fPartCross;             //![fMaxPerBasket] Indices array for crossing particles
   TGeoHelix        *fFieldPropagator;       //! Propagator in magnetic field
   TGeoRotation     *fRotation;              //! Rotation used by the field propagator
   GeantTrackCollection *fCollection;        //! Track collection per thread
   Int_t             fTracksPerBasket;

   std::vector<Int_t> fFinishedEvents;
	void PushFinishedEvent (Int_t eventSlot) { fFinishedEvents.push_back(eventSlot); }
	Int_t GetSizeOfFinishedEvents () const { return fFinishedEvents.size(); }

public:
   GeantThreadData();
   GeantThreadData(Int_t maxperbasket, Int_t maxprocesses);
   GeantThreadData(const GeantThreadData& oth);
   ~GeantThreadData();

   Double_t         *GetProcStep(Int_t iproc) {return fProcStep+iproc*fMaxPerBasket;}

};
#endif
