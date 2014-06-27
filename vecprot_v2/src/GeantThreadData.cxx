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
             fSizeDbl(0),
             fSizeBool(0),
             fVolume(0),
             fRndm(0),
             fBoolArray(0),
             fDblArray(0),
             fProcStep(0),
             fFieldPropagator(0),
             fRotation(0),
             fTrack(0),
             fPath(0)
{
// I/O ctor.
}

//______________________________________________________________________________
GeantThreadData::GeantThreadData(Int_t maxperbasket, Int_t maxprocesses)
            :TObject(),
             fMaxPerBasket(maxperbasket),
             fNprocesses(maxprocesses),
             fSizeDbl(5*maxperbasket),
             fSizeBool(5*maxperbasket),
             fVolume(0),
             fRndm(new TRandom()),
             fBoolArray(0),
             fDblArray(0),
             fProcStep(0),
             fFieldPropagator(0),
             fRotation(new TGeoRotation()),
             fTrack(0),
             fPath(0)
{
// Constructor
   GeantPropagator *propagator = GeantPropagator::Instance();
   fDblArray  = new Double_t[fSizeDbl];
   fBoolArray = new Bool_t[fSizeBool];
   fProcStep  = new Double_t[fNprocesses*fMaxPerBasket];
   fFieldPropagator = new TGeoHelix(1,1);
   fFieldPropagator->SetField(0,0,propagator->fBmag, kFALSE);
   fPath = new VolumePath_t(50); // hardcoded
}

//______________________________________________________________________________
GeantThreadData::~GeantThreadData()
{
// Destructor
 //  delete fMatrix;
   delete fRndm;
   delete [] fBoolArray;
   delete [] fDblArray;
   delete [] fProcStep;
   delete fFieldPropagator;
   delete fRotation;
   delete fPath;
}

//______________________________________________________________________________
Double_t *GeantThreadData::GetDblArray(Int_t size)
{
// Return memory storage for an array of doubles of at least "size" length which
// is thread-specific
   if (size<fSizeDbl) return fDblArray;
   Double_t *array = new Double_t[size];
   delete [] fDblArray;
   fDblArray = array;
   return fDblArray;
}

//______________________________________________________________________________
Bool_t *GeantThreadData::GetBoolArray(Int_t size)
{
// Return memory storage for an array of bool of at least "size" length which
// is thread-specific
   if (size<fSizeBool) return fBoolArray;
   Bool_t *array = new Bool_t[size];
   delete [] fBoolArray;
   fBoolArray = array;
   return fBoolArray;
}
