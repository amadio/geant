#ifndef GEANT_OUTPUT
#define GEANT_OUTPUT

#include "globals.h"

//______________________________________________________________________________
class GeantOutput : public TObject {
public:
   Double_t        fCpuTime;                 // Cpu time
   Int_t           fVolId;                   // Volume transporting this generation
   Int_t           fBasketGeneration;        // Burent generation of baskets to be flushed
   Int_t           fGeneration;              // Current generation for one basket
   Int_t           fNtracks;                 // Number of tracks in current generation
   Int_t          *fEvent;                   //[fNtracks]
   Int_t          *fInd;                     //[fNtracks] Track indices
   Int_t          *fProc;                    //[fNtracks] Selected processes for each track
   Double_t       *fX;                       //[fNtracks] X positions
   Double_t       *fY;                       //[fNtracks] Y positions
   Double_t       *fZ;                       //[fNtracks] Z positions
   Double_t       *fPx;                      //[fNtracks] Px
   Double_t       *fPy;                      //[fNtracks] Py
   Double_t       *fPz;                      //[fNtracks] Pz
   Double_t       *fE;                       //[fNtracks] E
   Double_t       *fPstep;                   //[fNtracks] Physics step selected
   Double_t       *fStep;                    //[fNtracks] Current step
   Double_t       *fSnext;                   //[fNtracks] Snext distance
   Double_t       *fSafety;                  //[fNtracks] Snext distance

private:
   GeantOutput(const GeantOutput&);  // Not implemented
   GeantOutput &operator=(const GeantOutput&);  // Not implemented

public:
   GeantOutput() : TObject(),fCpuTime(0),fVolId(-1),fBasketGeneration(0),fGeneration(0),fNtracks(0),fEvent(0),fInd(0),fProc(0),fX(0),fY(0),fZ(0),fPx(0),fPy(0),fPz(0),fE(0),fPstep(0),fStep(0),fSnext(0),fSafety(0) {}
   virtual ~GeantOutput();
   void            Init(Int_t size);
   void            Reset();
   void            SetStamp(Int_t volId, Int_t basket_gen, Int_t generation, Int_t ntracks, Double_t cputime=0.) {fVolId=volId; fBasketGeneration=basket_gen; fGeneration=generation;fNtracks=ntracks;fCpuTime=cputime;}
   void            SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t e, Double_t pstep, Double_t step, Double_t snext, Double_t safety);
   void            SetTrack(Int_t ntrack, GeantTrack *track);

   ClassDef(GeantOutput,1)       // The transport output per generation
};
#endif
