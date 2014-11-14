// Implementation - for source file

//  This class interfaces between the GeantV Propagator and 
//    the Compton model GUComptonKleinNishina
// 
// Authors:  J. Apostolakis, M. Novak 

#include "GUTrack.h"
#include "GUComptonKleinNishina.h"

#include "GVComptonProcess.h"

#include "globals.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"

FQUALIFIER
GVComptonProcess::GVComptonProcess(int processId) : 
  fProcessId (processId) 
{ 
  fProjectilesV= Create_GUTrack_v(MAX_NUMBER_TRACKS); 
  fSecondariesV= Create_GUTrack_v(MAX_NUMBER_TRACKS); 
  fpPhysicsModel= new GUComptonKleinNishina();
} 

FQUALIFIER
GVComptonProcess::~GVComptonProcess()
{ 
   delete fProjectilesV;
   delete fSecondariesV;
   delete fpPhysicsModel;
}

// This method should belong to a factory for GUTrack_v - or it is a class
GVComptonProcess::Create_GUTrack_v

FQUALIFIER int 
GVComptonProcess::ApplyPostStepProcess( GeantTrack_v& gTrackV )
{
  // Fill fProjectilesV // Must create it ...


  fpPhysicsModel->Interact( // int ntracks,
            fProjectilesV,    // In/Out
            fTargetElements,  // Number equal to num of tracks
            fSecondariesV);  // Empty vector for secondaries

  // Update status of primaries in gTrackV


  // Add secondaries to gTrackV - or where-ever

  for( i = ) 
  {	
     if(secEkin >= fEnergyLimit) //insert secondary into OUT tracks_v and rotate 
     {	
    	  GeantTrack &gTrack = propagator->GetTempTrack(tid);
//         GeantTrack gTrack;
        //set the new track properties
    	  gTrack.fEvent    = tracks.fEventV[t];
    	  gTrack.fEvslot   = tracks.fEvslotV[t];
        // gTrack.fParticle = nTotSecPart;          //index of this particle
	      gTrack.fPDG      = secPDG;               //PDG code of this particle
	      gTrack.fG5code   = pid[i];               //G5 index of this particle
	      gTrack.fEindex   = 0;
	      gTrack.fCharge   = secPartPDG->Charge()/3.; //charge of this particle
	      gTrack.fProcess  = 0;
	      gTrack.fIzero    = 0;
	      gTrack.fNsteps   = 0;
	//          gTrack.fSpecies  = 0;
	      gTrack.fStatus   = kNew;                 //status of this particle
	      gTrack.fMass     = secMass;              //mass of this particle
	      gTrack.fXpos     = tracks.fXposV[t];     //rx of this particle (same as parent)
	      gTrack.fYpos     = tracks.fYposV[t];     //ry of this particle (same as parent)
	      gTrack.fZpos     = tracks.fZposV[t];     //rz of this particle (same as parent)
	      gTrack.fXdir     = px/secPtot;     //dirx of this particle (before transform.)
	      gTrack.fYdir     = py/secPtot;     //diry of this particle before transform.)
	      gTrack.fZdir     = pz/secPtot;     //dirz of this particle before transform.)
	      gTrack.fP        = secPtot;              //momentum of this particle 
	      gTrack.fE        = secEtot;              //total E of this particle 
	      gTrack.fEdep     = 0.;
	      gTrack.fPstep    = 0.;
	      gTrack.fStep     = 0.;
	      gTrack.fSnext    = 0.;
	      gTrack.fSafety   = tracks.fSafetyV[t];
	      gTrack.fFrombdr  = tracks.fFrombdrV[t];
	      gTrack.fPending  = kFALSE;
	      *gTrack.fPath    = *tracks.fPathV[t];
	      *gTrack.fNextpath = *tracks.fPathV[t];

	      // Rotate new track to parent track's frame      
	      RotateNewTrack(oldXdir, oldYdir, oldZdir, gTrack);

	      propagator->AddTrack(gTrack);
	      tracks.AddTrack(gTrack);

	      ++nTotSecPart;
	  }    
}