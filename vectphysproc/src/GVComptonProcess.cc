// Implementation - for source file

//  This class interfaces between the GeantV Propagator and 
//    the Compton model GUComptonKleinNishina
// 
// Authors:  J. Apostolakis, M. Novak 
#include <iostream>

#include "GUTrack.h"
#include "GUTrackHandler.h"
#include "GUComptonKleinNishina.h"

#include "GVComptonProcess.h"

#include "globals.h"
#include "GeantTrack.h"
//#include "GeantPropagator.h"

#include "TPartIndex.h"

//ClassImp(GVComptonProcess)

FQUALIFIER
GVComptonProcess::GVComptonProcess(int processId, double energyLimit) : 
  fProcessId (processId), fEnergyLimit (energyLimit) 
{ 
  Create_GUTrack_v(&fProjectilesV,MAX_NUMBER_TRACKS); 
  Create_GUTrack_v(&fSecondariesV,MAX_NUMBER_TRACKS); 
  fpPhysicsModel= new GUComptonKleinNishina();
  std::cout<<"------------ process index: = "<< fProcessId << std::endl; 
} 


FQUALIFIER
GVComptonProcess::~GVComptonProcess()
{ 
   delete fProjectilesV;
   delete fSecondariesV;
   delete fpPhysicsModel;
}

// This method should belong to a factory for GUTrack_v - or it is a class
void GVComptonProcess::Create_GUTrack_v(GUTrack_v **tracks, int size){
  //GUTrackHandler *handler = new GUTrackHandler();
  //tracks = handler->GetSoATracks();
  *tracks = new GUTrack_v(size);
  
/*  dumy->numTracks = size;
    dumy->status    = new int(size);
  int *particleType;
  int *id;
  int *parentId;
  int *proc;
  double *x; 
  double *y;
  double *z;
  double *px;
  double *py;
  double *pz;
  double *E;
  double *q;
  double *s;
  double *fPdgV;        // PDG code
 */
 //*tracks = dumy;
}


