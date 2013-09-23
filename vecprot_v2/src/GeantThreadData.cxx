#include "GeantThreadData.h"
#include "globals.h"
#include "GeantBasket.h"
#include "GeantPropagator.h"

#include "TGeoMatrix.h"
#include "TArrayI.h"
#include "TGeoVolume.h"
#include "TRandom.h"
#include "TGeoHelix.h"

ClassImp(GeantThreadData)

//______________________________________________________________________________
GeantThreadData::GeantThreadData()
            :TObject(),
             fMaxPerBasket(0),
             fNprocesses(0),
             fMatrix(0),
             fVolume(0),
             fRndm(0),
             fDblArray(0),
             fProcStep(0),
             fFieldPropagator(0),
             fRotation(0)
{
// I/O ctor.
}

//______________________________________________________________________________
GeantThreadData::GeantThreadData(Int_t maxperbasket, Int_t maxprocesses)
            :TObject(),
             fMaxPerBasket(maxperbasket),
             fNprocesses(maxprocesses),
             fMatrix(new TGeoHMatrix()),
             fVolume(0),
             fRndm(new TRandom()),
             fDblArray(0),
             fProcStep(0),
             fFieldPropagator(0),
             fRotation(new TGeoRotation())
{
// Constructor
   GeantPropagator *propagator = GeantPropagator::Instance();
   fDblArray  = new Double_t[5*fMaxPerBasket];
   fProcStep  = new Double_t[fNprocesses*fMaxPerBasket];
   fFieldPropagator = new TGeoHelix(1,1);
   fFieldPropagator->SetField(0,0,propagator->fBmag, kFALSE);
}

//______________________________________________________________________________
GeantThreadData::~GeantThreadData()
{
// Destructor
   delete fMatrix;
   delete fRndm;
   delete [] fDblArray;
   delete [] fProcStep;
   delete fFieldPropagator;
   delete fRotation;
}
