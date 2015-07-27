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
   int             fMaxPerBasket;          //! Max number of tracks per basket
   int             fNprocesses;            //! Number of physics processes
   TGeoHMatrix      *fMatrix;                //! Current matrix per thread
   TGeoVolume       *fVolume;                //! Current volume per thread
   TRandom          *fRndm;                  //! Random generator for thread
   double         *fDblArray;              //![5*fMaxPerBasket] Support for random seed
   double         *fProcStep;              //![fNprocesses*fMaxPerBasket]
   TArrayI          *fPartInd;               //![fMaxPerBasket] Indices array for transported particles
   TArrayI          *fPartNext;              //![fMaxPerBasket] Indices array for remaining particles
   TArrayI          *fPartTodo;              //![fMaxPerBasket] Indices array for todo particles
   TArrayI          *fPartCross;             //![fMaxPerBasket] Indices array for crossing particles
   TGeoHelix        *fFieldPropagator;       //! Propagator in magnetic field
   TGeoRotation     *fRotation;              //! Rotation used by the field propagator
   GeantTrackCollection *fCollection;        //! Track collection per thread
   int             fTracksPerBasket;

   std::vector<int> fFinishedEvents;
	void PushFinishedEvent (int eventSlot) { fFinishedEvents.push_back(eventSlot); }
	int GetSizeOfFinishedEvents () const { return fFinishedEvents.size(); }

public:
   GeantThreadData();
   GeantThreadData(int maxperbasket, int maxprocesses);
   GeantThreadData(const GeantThreadData& oth);
   ~GeantThreadData();

   double         *GetProcStep(int iproc) {return fProcStep+iproc*fMaxPerBasket;}

};
#endif
