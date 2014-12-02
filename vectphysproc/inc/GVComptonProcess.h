#ifndef GVCOMPTONPROCESS_H
#define GVCOMPTONPROCESS_H

#define FQUALIFIER 
#define CONSTTYPE

#include "Geant/Config.h"

//#include<iostream>
// #include "GeantTrack.h"
//#include "GUTrack.h"
//#include "GUTrackHandler.h"
// #include "GUComptonKleinNishina.h"

class GeantTrack;
class GeantTrack_v;
class GUComptonKleinNishina;
class GUTrack_v;

class GVComptonProcess
{
  public:
  	FQUALIFIER GVComptonProcess();
  	FQUALIFIER GVComptonProcess(int processId, double energyLimit, 
                                    int maxnumtracks);
  	FQUALIFIER ~GVComptonProcess();
        // driver method for performing all the necessary actions
        // returns with the number of secondary tracks inserted into GeantTrack_v  
        FQUALIFIER int ApplyPostStepProcess(GeantTrack_v& gTrackV, int numtracks,
                                            int tid);
        
  private:
       void  Allocator(int size);
       void  Deallocator();
       // some GUTrack_v handling methods for allocation and de-allocation
       // these should be part of a GUTrack class later
       void  GUTrackAllocator(GUTrack_v& gutrack_v, int size);
       void  GUTrackDeallocator(GUTrack_v& gutrack_v);

       //  3-steps of the interaction:
       //# 3/1 -take the proper (i.e. Compton) primary tracks from the input GeantTrack_v 
       //      -set up the corresponding GUTrack_v         
       FQUALIFIER void FilterPrimaryTracks(GeantTrack_v& gTrackV , int numtracks);     
       //# 3/2 - call this vector physics model and perform the phyisics interaction
       FQUALIFIER void PerformInteraction();
       //# 3/3 - update primary tracks in GeantTrack_v to their post-interaction state
       //      - insert 'good' secondary tracks into GeantTrack_v 
       //      - return number of tracks inserted into GeantTrack_v
       FQUALIFIER int WriteBackTracks(GeantTrack_v& gTrackV, int tid);
       
       // initialize some members of a temporary GeantTrack 
       void SetGeantTrack(GeantTrack& left, GeantTrack_v& right, int ip);
                         
  private:
    // Parameters
    int     fProcessId;       // Index of this process in tabulated physics (TPartIndex)
                              // Used to select primary tracks for current interation from GeantTrack_v
    double  fEnergyLimit;     // Tracking cut in kinetic energy (In GeantV units [GeV]) 
    GUComptonKleinNishina*  fVComptonModel;  // this vector physics model 
     
    // Intermediate containers
    int     *fTargetElements;       // Z-of target atoms
    int     *fParentTrackIndices;   // indices of the parent tarcks in GeantTrack_v 
    GUTrack_v *fPrimaryTracks;      // GUTrack_v with the proper primary tracks for Compton 
    GUTrack_v *fSecondaryTracks;    // GUTrack_v with the corresponding secondary tracks 
                                    // will be filled by GUComptonKleinNishina::Interact

};  

#endif 

