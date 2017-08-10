#include "ExN03ApplicationRP.h"

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

#include "Geant/Error.h"

#include "PhysicsData.h"
#include "LightTrack.h"

#include <iostream>
#include <iomanip>

using std::min;
using std::max;

//______________________________________________________________________________
ExN03ApplicationRP::ExN03ApplicationRP(GeantRunManager *runmgr)
  : GeantVApplication(runmgr), fInitialized(false), fIdGap(0), fIdAbs(0), fNumWThreads(0), fFactory(0) {
  // Ctor..
  GeantFactoryStore *store = GeantFactoryStore::Instance();
  fFactory = store->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());
}

//______________________________________________________________________________
bool ExN03ApplicationRP::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  if (!GeoManager::Instance().GetWorld()) {
    Geant::Error("ExN03ApplicationRP::Initialize", "Geometry not loaded");
    return false;
  }
  Volume_t *lvGap = GeoManager::Instance().FindLogicalVolume("liquidArgon");
  Volume_t *lvAbs = GeoManager::Instance().FindLogicalVolume("Lead");
  if (!lvGap || !lvAbs) {
    Geant::Error("ExN03ApplicationRP::Initialize", "Logical volumes for gap and absorber not found - do you use the right geometry");
    return false;
  }
  fIdGap = lvGap->id();
  fIdAbs = lvAbs->id();
  //
  // set up the data structure to store data per working thread
  int numWThreads = fRunMgr->GetNthreadsTotal(); // number of working threads
  fNumWThreads    = numWThreads;
  fListDataPerThread.resize(numWThreads);
  int numThreads  = numWThreads + 5;                        // number of all threads (overestimate a bit)
  fWThreadIdToIndexMap.resize(numThreads);
  for (int i=0; i<numWThreads; ++i) {
    // first map the working thread ID to index
    int curWThreadId = fRunMgr->GetTDManager()->GetTaskData(i)->fTid;
    fWThreadIdToIndexMap[curWThreadId] = i;
    // set up the data structure for this thread with index i
    fListDataPerThread[i].fListDataPerAbsorber.resize(kNumAbsorbers);
    for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fEdep.resize(kNlayers,0.0);
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fLength.resize(kNlayers,0.0);
    }
  }

  fInitialized = true;
  return true;
}

//______________________________________________________________________________
void ExN03ApplicationRP::StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) {
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
    if (idvol==fIdAbs) { iabs = 0; }
    else if (idvol==fIdGap) { iabs = 1;}
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
      if (idvol==fIdGap && tracks.fEdepV[i]>0.00002) {
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


void ExN03ApplicationRP::SteppingActions(GeantTrack &track, GeantTaskData *td) {
  // Application stepping action.
  if (!fInitialized)
    return; // FOR NOW
  // energy deposit and step length
  int tid = td->fTid;
  Node_t const *current;
  int idvol = -1;
  int idnode = -1;
  int ilev = -1;
//  for (int i = 0; i < npart; i++) {
    ilev = track.Path()->GetCurrentLevel() - 1;
    if (ilev < 1)
      return;
    current = track.Path()->Top();
    if (!current)
      return;
    idnode = track.Path()->At(ilev - 1)->id();
    idvol  = current->GetLogicalVolume()->id();
    int indx   = fWThreadIdToIndexMap[tid];
    int ilayer = idnode;
    int iabs   = -1;
    if (idvol==fIdAbs) { iabs = 0; }
    else if (idvol==fIdGap) { iabs = 1;}
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
      if (idvol==fIdGap && track.fEdep>0.00002) {
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



void ExN03ApplicationRP::FinishRun() {
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
  int realNumLayesr = 10;
  int indxLayes[]   = {2,5,6,7,8,9,10,11,12,13};
  std::cout<< std::endl << std::setw(90) << std::setfill('=') << "" <<"\n";
  std::cout<< "  " << std::setw(30)  << std::setfill('*') << ""
           << "  RESULT OF THE SIMULATION "
           << std::setw(30) <<"" <<"\n";
  std::cout<< std::setw(90) << std::setfill('=') << "" <<std::setfill(' ') <<"\n";
  std::cout<<" [----- Energy deposit [MeV/#primary] per layer: [#Layer]  [Lead]  [Scintillator] ------] \n";
  for (int i=0; i<realNumLayesr; ++i) {
    int ilayer = indxLayes[i];
    std::cout<< std::right << std::setw(4) << i
             << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[0].fEdep[ilayer]*1000./nprim
             << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[1].fEdep[ilayer]*1000./nprim
             << std::endl;
  }
  std::cout << std::endl;
  std::cout<<" [-- Cumulated track length [cm/#primary] per layer: [#Layer]  [Lead]  [Scintillator] --] \n";
  for (int i=0; i<realNumLayesr; ++i) {
    int ilayer = indxLayes[i];
    std::cout<< std::right << std::setw(4) << i
             << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[0].fLength[ilayer]/nprim
             << std::setw(13) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[1].fLength[ilayer]/nprim
             << std::endl;
  }
  std::cout<< std::setw(90) << std::setfill('=') << "" << "\n\n";
}
