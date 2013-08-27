#include "GeantHitBlock.h"
#include "TClass.h"

ClassImp(GeantHitBlock)

//______________________________________________________________________________
GeantHitBlock::GeantHitBlock(TClass *cl,Int_t size, void *storage)
          :TObject(),
           fSize(size),
           fNext(0),
           fClass(cl),
           fClSize(0),
           fSupport(0),
           fArray(0)
{
// Constructor
   if (!cl) throw "Null class provided";
   fClSize = cl->Size();
   fSupport = ( storage ) ? new (storage) char[size*fClSize+8] : new char[size*fClSize+8];
   if (!fSupport) throw "Allocation error for GeantHitBlock";
   fArray = cl->NewArray(size, fSupport);
}

//______________________________________________________________________________
GeantHitBlock::~GeantHitBlock()
{
// Destructor
  for (Int_t i=0; i<fSize; i++) {
     TObject *obj = reinterpret_cast<TObject *>((ULong_t)fArray+fClSize*i);
     TObject::SetDtorOnly(obj);
     fClass->Destructor(obj);
  }
  // Delete the support array
  delete [] fSupport;
}

//______________________________________________________________________________
TObject* GeantHitBlock::NextHit()
{
// Get next free hit location. Returns 0 if the hit block is full
   if (fNext==fSize) return 0;
   return reinterpret_cast<TObject *>((ULong_t)fArray+fClSize*(fNext++));
}

//______________________________________________________________________________
void GeantHitBlock::ClearHits()
{
// Clear the hits, provided that the hit class overwrites 
// TObject::Clear(Option_t *option)
   for (Int_t i=0; i<fSize; i++) {
      TObject *obj = reinterpret_cast<TObject *>((ULong_t)fArray+fClSize*i);
      obj->Clear();
   }   
}

//______________________________________________________________________________
void GeantHitBlock::CopyHit(Int_t ihit, GeantHitBlock &other)
{
// Copy the hit into a different hit block. No check made (the other hit block
// should have free hits)
   TObject *hit = reinterpret_cast<TObject *>((ULong_t)fArray+fClSize*ihit);
   TObject *dest = other.NextHit();
      
}
