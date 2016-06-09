#include "globals.h"
#include "GeantOutput.h"
#include "GeantTrack.h"

//______________________________________________________________________________
GeantOutput::~GeantOutput() {
  // Destructor
  Reset();
}

//______________________________________________________________________________
void GeantOutput::Init(int size) {
  // Initialize arrays to a given size.
  Reset();
  fEvent = new int[size];
  fInd = new int[size];
  fProc = new int[size];
  fX = new double[size];
  fY = new double[size];
  fZ = new double[size];
  fPx = new double[size];
  fPy = new double[size];
  fPz = new double[size];
  fE = new double[size];
  fPstep = new double[size];
  fStep = new double[size];
  fSnext = new double[size];
  fSafety = new double[size];
}

//______________________________________________________________________________
void GeantOutput::Reset() {
  // Reset arrays
  delete[] fEvent;
  fEvent = 0;
  delete[] fInd;
  fInd = 0;
  delete[] fProc;
  delete[] fX;
  delete[] fY;
  delete[] fZ;
  delete[] fPx;
  delete[] fPy;
  delete[] fPz;
  delete[] fE;
  delete[] fPstep;
  delete[] fStep;
  delete[] fSnext;
  delete[] fSafety;
}

//______________________________________________________________________________
void GeantOutput::SetTrack(int ntrack, int itrack, int event, int proc, double x,
                           double y, double z, double px, double py, double pz,
                           double e, double pstep, double step, double snext,
                           double safety) {
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
void GeantOutput::SetTrack(int ntrack, GeantTrack *track) {
  // Set parameters for ntrack based on a GeantTrack
  SetTrack(ntrack, track->fParticle, track->fEvent, track->fProcess, track->fXpos, track->fYpos,
           track->fZpos, track->Px(), track->Py(), track->Pz(), track->fE, track->fPstep,
           track->fStep, track->fSnext, track->fSafety);
}
