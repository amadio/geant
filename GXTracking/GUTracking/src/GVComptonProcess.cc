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
  int numTracksIn = gTrackV.number;
  int num= 0;

  for ( int i=0; i<numTracksIn; i++)   // Assumes that there are no holes
  {
  	 if( gTrackV.fProcessV[i] == fProcessId)
  	 { 
  	 	num++;
  	 	fIndex[num]=  i;
        // fProjectilesV.status[num]=  gTrackV.fStatusV[i];   // Care: type is TrackStatus_t 
        // fProjectilesV.particleType[num]=  gTrackV.fParticleV / fG5codeV / fPDGV / Species_t *fSpeciesV;
        fProjectilesV.id[num]=       num;  // gTrackV.fParticleV[i];  //  or num - or i ?
        // fProjectilesV.parentId[num]= gTrackV.[i];        // Not available - ok
        // fProjectilesV.proc[num]=  gTrackV.fProcessV[i];  // Process Index -- this one 

        fProjectilesV.x[num]=  gTrackV.fXposV[i];      // Just to check it - not really needed
        fProjectilesV.y[num]=  gTrackV.fYposV[i];
        fProjectilesV.z[num]=  gTrackV.fZposV[i];

        double momentum=  gTrackV.fPV[i];
        fProjectilesV.px[num]=  momentum*gTrackV.fXdirV[i];
        fProjectilesV.py[num]=  momentum*gTrackV.fYdirV[i];
        fProjectilesV.pz[num]=  momentum*gTrackV.fZdirV[i];
        fProjectilesV.E[num]=  gTrackV.fEV[i];
        // fProjectilesV.q[num]=  gTrackV.fChargeV[i];
        // fProjectilesV.s[num]=  gTrackV.fPStepV[i];

        fTargetElements[num]= fEindexV[num];  // Mihaly must check this is Z -- else fix TODO
        // =gTrackV.fEindexV[i];
     }
  }
  fProjectilesV.numTracks= num;

  fpPhysicsModel->Interact( // int ntracks,
            fProjectilesV,    // In/Out
            fTargetElements,  // Z of elements. [ Size equal to num of tracks ]
            fSecondariesV);   // Empty vector for secondaries

  // Update status of primaries in gTrackV using fProjectilesV
  for ( int j=0; j<num; j++)
  { 
  	int i= fIndex[j];

    G4double  Etot= fProjectilesV.E[j];
    G4double  Ekin= Etot- gTrackV.fMassV[i];  // Must be in GeV
    gTrackV.fEV[i]= Etot;     // fProjectilesV.E[j];
    if( Ekin < fEnergyLimit) {
    	gTrackV.fEdepV[i] += Ekin;
    	gTrackV.fStatusV[i] = Ekilled;  // Fix - ToDo MN 	
    }else{
        // gTrackV.[i]= fProjectilesV.parentId[i]= ;
        // fProjectilesV.proc[i]=  gTrackV.[i];
        // gTrackV.fXposV[i] = fProjectilesV.x[i];
        // gTrackV.fYposV[i] = fProjectilesV.y[i];
        // gTrackV.fZposV[i] = fProjectilesV.z[i];
        G4double px= fProjectilesV.px[j];
        G4double py= fProjectilesV.py[j];
        G4double pz= fProjectilesV.pz[j];
        // G4double momentumSq= (px*px+py*py+pz*pz);
        // G4double momentum= std::sqrt(momentumSq);
        // G4double invMomentum= 1.0 / momentum;
        G4double invMomentum= 1.0 / std::sqrt(px*px+py*py+pz*pz);
        gTrackV.fXdirV[i]= invMomentum * px;
        gTrackV.fYdirV[i]= invMomentum * py;
        gTrackV.fZdirV[i]= invMomentum * pz;
    }
    // fProjectilesV.q[i]
  }

  // Add secondaries to gTrackV - or where-ever
  int numCandSecondaries= fSecondariesV.numTracks;

  int numInserted = 0;   // Counter for secondaries that survive production cut

  const TDatabasePDG *dbPdg= TPartIndex::I()->DBPdg();

  // A process with the same secondary can define its base properties here
  const int          secPDG =        22 ;    // PDG code of photon
  const int          secG5code=      TPartIndex::I()->PartIndex(secPDG); 
  const TParticlePDG secParticleDef= dbPdb->...
  const double       secMass = 0.0;  // photon rest mass
  const double       secCharge = 0; 

  for ( int iSec=0; iSec<numCandSecondaries; iSec++)
  {	
  	 int    idParent = fSecondariesV.parentId[iSec];
     // fIndex[iSec];
     double secEtot= fSecondariesV.fEV[i];
     // const int    secPDG = fSecondariesV.fPdgV[i];
     // const int    secG5code= TPartIndex::I()->PartIndex(secPDG); 
     // const double secMass = secParticleDef->GetMass() ....   // TODO

     if(secEkin < fEnergyLimit) {
        gTrackV.fEdepV[idParent] += secEtot - secMass ;
     }
     else //insert secondary into OUT tracks_v and rotate 
     {	
          double px, py, pz; 
          px= fSecondariesV.px[iSec];
          py= fSecondariesV.py[iSec];
          pz= fSecondariesV.pz[iSec];                   
          double secPtot= std::sqrt(px*px+py*py+pz*pz); 
          double inv_Ptot= 1.0 / secPtot;

          // gTrackV.fParticleV[i] fProjectilesV.id[num]= ;  //  or num - or i ?     	
    	  GeantTrack &gTrack = propagator->GetTempTrack(tid);
          // GeantTrack gTrack;
        //set the new track properties
    	  gTrack.fEvent    = fProjectilesV.fEventV[idParent];
    	  gTrack.fEvslot   = fProjectilesV.fEvslotV[idParent];
        // gTrack.fParticle = nTotSecPart;          //index of this particle
	      gTrack.fPDG      = secPDG;               //PDG code of this particle
	      gTrack.fG5code   = TPartIndex::I()->PartIndex(secPDG);  //G5 index of this particle
	      gTrack.fEindex   = 0;
	      gTrack.fCharge   = secCharge; //charge of this particle
	      gTrack.fProcess  = 0;
	      gTrack.fIzero    = 0;
	      gTrack.fNsteps   = 0;
	   // gTrack.fSpecies  = 0;
	      gTrack.fStatus   = kNew;                 //status of this particle
	      gTrack.fMass     = secMass;              //mass of this particle
	      gTrack.fXpos     = gTrackV.fXposV[idParent];     //rx of this particle (same as parent)
	      gTrack.fYpos     = gTrackV.fYposV[idParent];     //ry of this particle (same as parent)
	      gTrack.fZpos     = gTrackV.fZposV[idParent];     //rz of this particle (same as parent)
	      gTrack.fXdir     = px*inv_Ptot;     //dirx of this particle (before transform.)
	      gTrack.fYdir     = py*inv_Ptot;     //diry of this particle before transform.)
	      gTrack.fZdir     = pz*inv_Ptot;     //dirz of this particle before transform.)
	      gTrack.fP        = secPtot;              //momentum of this particle 
	      gTrack.fE        = secEtot;              //total E of this particle 
	      gTrack.fEdep     = 0.;
	      gTrack.fPstep    = 0.;
	      gTrack.fStep     = 0.;
	      gTrack.fSnext    = 0.;
	      gTrack.fSafety   = gTrackV.fSafetyV[idParent];
	      gTrack.fFrombdr  = gTrackV.fFrombdrV[idParent];
	      gTrack.fPending  = kFALSE;
	      *gTrack.fPath    = *gTrackV.fPathV[idParent];
	      *gTrack.fNextpath = *gTrackV.fPathV[idParent];

	      propagator->AddTrack(gTrack);
	      gTrackV.AddTrack(gTrack);
          ++numInserted;
	  }
      return numInserted; 
}