#include "globals.h"
#include "GeantOutput.h"
#include "GeantTrack.h"

//ClassImp(GeantOutput)

//______________________________________________________________________________
GeantOutput::~GeantOutput()
{
// Destructor
   Reset();
}

//______________________________________________________________________________
void GeantOutput::Init(Int_t size)
{
// Initialize arrays to a given size.
   Reset();
   fEvent = new Int_t[size];
   fInd = new Int_t[size];
   fProc = new Int_t[size];
   fX = new Double_t[size];
   fY = new Double_t[size];
   fZ = new Double_t[size];
   fPx = new Double_t[size];
   fPy = new Double_t[size];
   fPz = new Double_t[size];
   fE  = new Double_t[size];
   fPstep = new Double_t[size];
   fStep = new Double_t[size];
   fSnext = new Double_t[size];
   fSafety = new Double_t[size];
}   
   
//______________________________________________________________________________
void GeantOutput::Reset()
{
// Reset arrays
   delete [] fEvent; fEvent = 0;
   delete [] fInd; fInd = 0;
   delete [] fProc;
   delete [] fX; delete [] fY; delete [] fZ;
   delete [] fPx; delete [] fPy; delete [] fPz; delete [] fE;
   delete [] fPstep; delete [] fStep; delete [] fSnext; delete [] fSafety;
}   

//______________________________________________________________________________
void GeantOutput::SetTrack(Int_t ntrack, Int_t itrack, Int_t event, Int_t proc, Double_t x, Double_t y, Double_t z, Double_t px, Double_t py, Double_t pz, Double_t e, Double_t pstep, Double_t step, Double_t snext, Double_t safety)
{
// Set parameters for ntrack
   fInd[ntrack] = itrack;
   fEvent[ntrack] = event;
   fProc[ntrack] = proc;
   fX[ntrack] = x;
   fY[ntrack] = y;
   fZ[ntrack] = z;
   fPx[ntrack] = px;
   fPy[ntrack] = py;
   fPz[ntrack] = pz;
   fE[ntrack] = e;
   fPstep[ntrack] = pstep;
   fStep[ntrack] = step;
   fSnext[ntrack] = snext;
   fSafety[ntrack] = safety;
}   

//______________________________________________________________________________
void GeantOutput::SetTrack(Int_t ntrack, GeantTrack *track)
{
// Set parameters for ntrack based on a GeantTrack
   SetTrack(ntrack, track->particle, track->event, track->process, track->xpos, track->ypos, track->zpos, track->px, track->py, track->pz, track->e, track->pstep, track->step, track->snext, track->safety);
}   
