///////////////////////////////////////////////////////////////////////////////////////////
////
////                      CaloApplicationRP.cxx
////                      Created: 20 June 2017
////                      Author: Ryan Schmitz
////
//// Description: A (linear) calorimeter implemented using VecGeom libraries. This construction is 
//			very similar to CaloApp.cxx, but is implemented using a different data structure
//			(outputting energy/track legnth per absorber, rather than per absorber type).
//			This example also doesn't collect true standard deviations; the newer CaloApp
//			example must be used to collect this information.
////
////////////////////////////////////////////////////////////////////////////////////////////


#include "CaloApplicationRP.h"

//#ifdef USE_VECGEOM_NAVIGATOR
#include "management/GeoManager.h"
using vecgeom::GeoManager;
//#endif
#include "GeantEvent.h"
#include "GeantFactoryStore.h"
#include "GeantTrackVec.h"
#include "GeantRunManager.h"
#include "GeantTaskData.h"
#include "globals.h"
#ifdef USE_ROOT
#include "TGeoNode.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#endif
#include <cassert>
#include "CaloDetectorConstruction.h"
#include "Geant/Error.h"

#include "PhysicsData.h"
#include "LightTrack.h"
#include <vector>
#include <iostream>
#include <iomanip>

using std::min;
using std::max;

//______________________________________________________________________________
CaloApplicationRP::CaloApplicationRP(GeantRunManager *runmgr,CaloDetectorConstruction *userCalo)
  : GeantVApplication(runmgr), fInitialized(false), fNumWThreads(0), fFactory(0) {
  // Ctor..
  GeantFactoryStore *store = GeantFactoryStore::Instance();
  fFactory = store->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());
  calo=userCalo;
}

//______________________________________________________________________________
bool CaloApplicationRP::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  if (!GeoManager::Instance().GetWorld()) {
    Geant::Error("CaloApplicationRP::Initialize", "Geometry not loaded");
    return false;
  }
//  vecgeom::cxx::LogicalVolume *lvGap[kNlayers];
  kNlayers = calo->GetNumLayers();
  kNumAbsorbers = calo->GetNumAbsorbers();
  std::cout << "Does this work? num abs = " << kNumAbsorbers << " and num layers = " << kNlayers << std::endl;
  vecgeom::cxx::LogicalVolume *lvAbs[kNlayers+1][kNumAbsorbers+1];
//  char *gapName = new char;
  char *absName = new char;
  for (int i=1;i<=kNlayers;i++){
	for(int j=1;j<=kNumAbsorbers; j++){
//  	sprintf(gapName,"abs%d1",i);
  	sprintf(absName,"abs%d%d",i,j);
//  	lvGap[i] = GeoManager::Instance().FindLogicalVolume(gapName);
  	lvAbs[i][j] = GeoManager::Instance().FindLogicalVolume(absName);
  	if (!lvAbs[i][j]) {
    		Geant::Error("CaloApplicationRP::Initialize", "Logical volumes for absorber not found - did you use the right geometry?");
   		return false;
  	}
//  fIdGap[i-1] = lvGap[i]->id();
  fIdAbs[i-1][j-1] = lvAbs[i][j]->id();
  }	
  }
  
  //
  // set up the data structure to store data per working thread
  int numWThreads = fRunMgr->GetNthreadsTotal(); // number of working threads
  fNumWThreads    = numWThreads;
  fListDataPerThread.resize(numWThreads);
  int numThreads  = numWThreads + 5;                        // number of all threads (overestimate a bit)
  fWThreadIdToIndexMap.resize(numThreads);
  for (int i=0; i<numWThreads; ++i) {
    // first map the working thread ID to index
    int curWThreadId = fRunMgr->GetTaskData(i)->fTid;
    fWThreadIdToIndexMap[curWThreadId] = i;
    // set up the data structure for this thread with index i
    fListDataPerThread[i].fListDataPerAbsorber.resize(kNumAbsorbers);
    for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fEdep.resize(kNlayers,0.0);
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fLength.resize(kNlayers,0.0);
    }
  }


//
// CREATE PhysicsData here: should be done at the init of PhysicsProcessHandler but
// GeantTaskData are constructed later than that call
  for (int i=0; i<fRunMgr->GetNthreadsTotal(); ++i) {
    fRunMgr->GetTaskData(i)->fPhysicsData = new geantphysics::PhysicsData();
  }

  fInitialized = true;
  return true;
};
//______________________________________________________________________________
void CaloApplicationRP::StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  if (!fInitialized)
    return; // FOR NOW
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  int tid = td->fTid;
  Node_t const *current;
  int idvol = -1;
  int idnode = -1;
  int ilev = -1;
  bool scintHit=false;
  for (int i = 0; i < npart; i++) {
    ilev = tracks.fPathV[i]->GetCurrentLevel() - 1;
    if (ilev < 1)
      continue;
    current = tracks.fPathV[i]->Top();
    if (!current)
      continue;
    idnode = tracks.fPathV[i]->At(ilev - 1)->id();
    idvol  = current->GetLogicalVolume()->id();
    int indx   = fWThreadIdToIndexMap[tid];
    int ilayer = idnode;
    int iabs   = -1;
    for (int j=0;j<kNlayers;j++){
        for (int k=0; k<kNumAbsorbers; k++){
    		if (idvol==fIdAbs[j][k]) {
			ilayer=j;
			iabs = k;
			if (iabs==0) scintHit=true;
			break;
                }
        }
    } 
    if (iabs>-1) {
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fEdep[ilayer]   +=  tracks.fEdepV[i];
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fLength[ilayer] +=  tracks.fStepV[i];
    }
  }

  if (fRunMgr->GetConfig()->fFillTree) {
    MyHit *hit;
    //    int nhits = 0;
    for (int i = 0; i < npart; i++) {
      // Deposit hits in scintillator
      if (scintHit && tracks.fEdepV[i]>0.00002) {
	      hit = fFactory->NextFree(tracks.fEvslotV[i], tid);
	      hit->fX = tracks.fXposV[i];
	      hit->fY = tracks.fYposV[i];
	      hit->fZ = tracks.fZposV[i];
	      hit->fEdep = tracks.fEdepV[i];
        hit->fTime = tracks.fTimeV[i];
        hit->fEvent = tracks.fEventV[i];
        hit->fTrack = tracks.fParticleV[i];
	      hit->fVolId = idvol;
	      hit->fDetId = idnode;
    	  //      if (track->path && track->path->GetCurrentNode()) {
    	  //         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
    	  //         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
    	  //      }
    	  //	  nhits++;
     	}
    }
  }

  //  Printf("Thread %d produced %d hits", tid, nhits);
  //  Printf("StepManager: size of queue %zu", fFactory->fOutputs.size());

  return;
}


void CaloApplicationRP::SteppingActions(GeantTrack &track, GeantTaskData *td) {
  // Application stepping action.
  if (!fInitialized)
    return; // FOR NOW
  // energy deposit and step length
  int tid = td->fTid;
  Node_t const *current;
  int idvol = -1;
  int idnode = -1;
  int ilev = -1;
  bool scintHit=false;
//  for (int i = 0; i < npart; i++) {
    ilev = track.fPath->GetCurrentLevel() - 1;
    if (ilev < 1)
      return;
    current = track.fPath->Top();
    if (!current)
      return;
    idnode = track.fPath->At(ilev - 1)->id();
    idvol  = current->GetLogicalVolume()->id();
    int indx   = fWThreadIdToIndexMap[tid];
    int ilayer = idnode;
    int iabs   = -1;
    for (int j=0;j<kNlayers;j++){
        for (int k=0; k<kNumAbsorbers; k++){
    		if (idvol==fIdAbs[j][k]) {
			ilayer=j;
			iabs = k;
			if (iabs==0) scintHit=true;
			break;
                }
        }
    } 
    if (iabs>-1) {
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fEdep[ilayer]   +=  track.fEdep;
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fLength[ilayer] +=  track.fStep;
    }
//  }


  if (fRunMgr->GetConfig()->fFillTree) {
    MyHit *hit;
    //    int nhits = 0;
//    for (int i = 0; i < npart; i++) {
      // Deposit hits in scintillator
      if (scintHit && track.fEdep>0.00002) {
	      hit = fFactory->NextFree(track.fEvslot, tid);
	      hit->fX = track.fXpos;
	      hit->fY = track.fYpos;
	      hit->fZ = track.fZpos;
	      hit->fEdep = track.fEdep;
        hit->fTime = track.fTime;
        hit->fEvent = track.fEvent;
        hit->fTrack = track.fParticle;
	      hit->fVolId = idvol;
	      hit->fDetId = idnode;
    	  //      if (track->path && track->path->GetCurrentNode()) {
    	  //         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
    	  //         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
    	  //      }
    	  //	  nhits++;
     	}
//    }
  }
}



void CaloApplicationRP::FinishRun() {
  double nprim = (double)fRunMgr->GetNprimaries();
  for (int i=0; i<kNlayers; ++i) {
    for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
      for (int iwthrd=1; iwthrd<fNumWThreads; ++iwthrd) {
        fListDataPerThread[0].fListDataPerAbsorber[iabs].fEdep[i]
        += fListDataPerThread[iwthrd].fListDataPerAbsorber[iabs].fEdep[i];
        fListDataPerThread[0].fListDataPerAbsorber[iabs].fLength[i]
        +=fListDataPerThread[iwthrd].fListDataPerAbsorber[iabs].fLength[i];
      }
    }
  }
//  int realNumLayesr = 10;
//  int indxLayes[]   = {0,1,2,3,4,5,6,7,8,9};
  std::cout<< std::endl << std::setw(90) << std::setfill('=') << "" <<"\n";
  std::cout<< "  " << std::setw(30)  << std::setfill('*') << ""
           << "  RESULT OF THE SIMULATION "
           << std::setw(30) <<"" <<"\n";
  std::cout<< std::setw(90) << std::setfill('=') << "" <<std::setfill(' ') <<"\n";
  std::cout<<" [--- Energy deposit [MeV/#primary] per layer ---]\n[#Layer]";
  for (int p=1; p<=kNumAbsorbers; p++)  std::cout<< "["<<calo->GetAbsorberMaterialName(p) << "]";
  std::cout<<"\n";
  for (int i=0; i<kNlayers; ++i) {
//    int ilayer = indxLayes[i];
    std::cout<< std::right << std::setw(4) << i;
	for (int j=0; j<kNumAbsorbers; j++){
             std::cout << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[j].fEdep[i]*1000./nprim;
             }
		std::cout << std::endl;
  }
  std::cout << std::endl;
  std::cout<<"[-- Cumulated track length [cm/#primary] per layer --]\n[#Layer]";
  for (int p=1; p<=kNumAbsorbers; p++)  std::cout<< "["<<calo->GetAbsorberMaterialName(p) << "]";
  std::cout<<"\n";
  for (int i=0; i<kNlayers; ++i) {
//    int ilayer = indxLayes[i];
    std::cout<< std::right << std::setw(4) << i;
	for (int j=0; j<kNumAbsorbers; j++){    
             std::cout << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[j].fLength[i]/nprim;
	}
             std::cout << std::endl;
  }
  std::cout<< std::setw(90) << std::setfill('=') << "" << "\n\n";
}
