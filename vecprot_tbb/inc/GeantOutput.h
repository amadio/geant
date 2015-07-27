#ifndef GEANT_OUTPUT
#define GEANT_OUTPUT

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class GeantTrack;

//______________________________________________________________________________
class GeantOutput : public TObject {
public:
   double        fCpuTime;                 // Cpu time
   Int_t           fVolId;                   // Volume transporting this generation
   Int_t           fGeneration;              // Current generation for one basket
   Int_t           fNtracks;                 // Number of tracks in current generation
   Int_t          *fEvent;                   //[fNtracks]
   Int_t          *fInd;                     //[fNtracks] Track indices
   Int_t          *fProc;                    //[fNtracks] Selected processes for each track
   double       *fX;                       //[fNtracks] X positions
   double       *fY;                       //[fNtracks] Y positions
   double       *fZ;                       //[fNtracks] Z positions
   double       *fPx;                      //[fNtracks] Px
   double       *fPy;                      //[fNtracks] Py
   double       *fPz;                      //[fNtracks] Pz
   double       *fE;                       //[fNtracks] E
   double       *fPstep;                   //[fNtracks] Physics step selected
   double       *fStep;                    //[fNtracks] Current step
   double       *fSnext;                   //[fNtracks] Snext distance
   double       *fSafety;                  //[fNtracks] Snext distance

public:
   GeantOutput() : TObject(),fCpuTime(0),fVolId(-1),fGeneration(0),
                     fNtracks(0),fEvent(0),fInd(0),fProc(0),fX(0),fY(0),fZ(0),fPx(0),fPy(0),fPz(0),
                     fE(0),fPstep(0),fStep(0),fSnext(0),fSafety(0) {}
   virtual ~GeantOutput();
   void            Init(Int_t size);
   void            Reset();
   void            SetStamp(Int_t volId, Int_t basket_gen, Int_t generation, Int_t ntracks, double cputime=0.)
                      {fVolId=volId; fGeneration=generation;fNtracks=ntracks;fCpuTime=cputime;}
   void            SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc,
                              double x, double y, double z, double px, double py, double pz,
                              double e, double pstep, double step, double snext, double safety);
   void            SetTrack(Int_t ntrack, GeantTrack *track);

   ClassDef(GeantOutput,1)       // The transport output per generation
};
#endif
