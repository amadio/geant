#ifndef GVCOMPTONPROCESS_H
#define GVCOMPTONPROCESS_H

#define FQUALIFIER 
#define CONSTTYPE

#include "Geant/Config.h"

#include<iostream>
// #include "GeantTrack.h"
#include "GUTrack.h"
//#include "GUTrackHandler.h"
// #include "GUComptonKleinNishina.h"

class GeantTrack;
class GeantTrack_v;
class GUComptonKleinNishina;
class GUTrack;
class GUTrackHandler;

#define  MAX_NUMBER_TRACKS   256

class GVComptonProcess
{
  public:
  	FQUALIFIER GVComptonProcess(int processId, double energyLimit);
  	FQUALIFIER ~GVComptonProcess();
 
      //  FQUALIFIER int ApplyPostStepProcess( GeantTrack_v& gTrackV ); // In/Out vector of tracks to apply
                                        
      // Returns the number of secondaries created

  protected:
    FQUALIFIER void Create_GUTrack_v(GUTrack_v **tracks, int size);

  private:
    // Parameters
    int     fProcessId;       // Used to select primaries for current interation
    double  fEnergyLimit;     // In GeantV units - now GeV 
    GUComptonKleinNishina* fpPhysicsModel;
    int     fNumPrimaries;    // 
    int     fNumSecondaries;  // Number created by our process - valid after call
     
  	// Intermediate containers - Output secondaries
    int         fTargetElements[MAX_NUMBER_TRACKS];
    GUTrack_v *fProjectilesV;
    GUTrack_v *fSecondariesV;
//    GUTrackHandler *fTrackHandelIn;
//    GUTrackHandler *fTrackHandelOut;

// ClassDef(GVComptonProcess,1) //
};  

#endif 

