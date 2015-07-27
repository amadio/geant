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
   int           fVolId;                   // Volume transporting this generation
   int           fGeneration;              // Current generation for one basket
   int           fNtracks;                 // Number of tracks in current generation
   int          *fEvent;                   //[fNtracks]
   int          *fInd;                     //[fNtracks] Track indices
   int          *fProc;                    //[fNtracks] Selected processes for each track
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
   void            Init(int size);
   void            Reset();
   void            SetStamp(int volId, int basket_gen, int generation, int ntracks, double cputime=0.)
                      {fVolId=volId; fGeneration=generation;fNtracks=ntracks;fCpuTime=cputime;}
   void            SetTrack(int ntrack, int itrack, int event, int proc,
                              double x, double y, double z, double px, double py, double pz,
                              double e, double pstep, double step, double snext, double safety);
   void            SetTrack(int ntrack, GeantTrack *track);

   ClassDef(GeantOutput,1)       // The transport output per generation
};
#endif
