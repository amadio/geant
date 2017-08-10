#include "TestEM3ApplicationRP.h"

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
#include "Particle.h"

#include <iostream>
#include <iomanip>

using std::min;
using std::max;

//______________________________________________________________________________
TestEM3ApplicationRP::TestEM3ApplicationRP(GeantRunManager *runmgr)
  : GeantVApplication(runmgr), fInitialized(false), fIdGap(0), fIdAbs(0), fNumWThreads(0), fFactory(0) {
  // Ctor..
  GeantFactoryStore *store = GeantFactoryStore::Instance();
  fFactory = store->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());
}

//______________________________________________________________________________
bool TestEM3ApplicationRP::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  if (!GeoManager::Instance().GetWorld()) {
    Geant::Error("TestEM3ApplicationRP::Initialize", "Geometry not loaded");
    return false;
  }
  Volume_t *lvGap = GeoManager::Instance().FindLogicalVolume("liquidArgon");
  Volume_t *lvAbs = GeoManager::Instance().FindLogicalVolume("Lead");
  if (!lvGap || !lvAbs) {
    Geant::Error("TestEM3ApplicationRP::Initialize", "Logical volumes for gap and absorber not found - do you use the right geometry");
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
    // 1. init global data i.e. not per absorber data
    fListDataPerThread[i].fDataGlobal.fNumGamma        = 0.0;
    fListDataPerThread[i].fDataGlobal.fNumElectron     = 0.0;
    fListDataPerThread[i].fDataGlobal.fNumPositron     = 0.0;
    fListDataPerThread[i].fDataGlobal.fNumChargedSteps = 0.0;
    fListDataPerThread[i].fDataGlobal.fNumNeutralSteps = 0.0;
    // 2. per absorber data
    fListDataPerThread[i].fListDataPerAbsorber.resize(kNumAbsorbers);
    for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fEdep   = 0.0;
      fListDataPerThread[i].fListDataPerAbsorber[iabs].fLength = 0.0;
    }
  }

  fInitialized = true;
  return true;
}


//______________________________________________________________________________
void TestEM3ApplicationRP::StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) {
  // Application stepping action.
  if (!fInitialized)
    return; // FOR NOW
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
//    int ilayer = idnode;
    int iabs   = -1;
    if (idvol==fIdAbs) { iabs = 0; }
    else if (idvol==fIdGap) { iabs = 1;}

    const geantphysics::Particle *part = geantphysics::Particle::GetParticleByInternalCode(tracks.fGVcodeV[i]);
    int   pdgCode = part->GetPDGCode();
    double charge = part->GetPDGCharge();

    if (iabs>-1) {
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fEdep     +=  tracks.fEdepV[i];
      if (charge!=0.0) {
        fListDataPerThread[indx].fListDataPerAbsorber[iabs].fLength +=  tracks.fStepV[i];
      }
    }
    if (tracks.fStatusV[i]!=Geant::kNew) { // do not count the creation step
      if (charge==0.0) {
        fListDataPerThread[indx].fDataGlobal.fNumNeutralSteps +=1.0;
      } else {
        fListDataPerThread[indx].fDataGlobal.fNumChargedSteps +=1.0;
      }
    }
    if (tracks.fStatusV[i]==Geant::kNew) { // check secondaries
      switch(pdgCode) {
        // e+
        case -11 : fListDataPerThread[indx].fDataGlobal.fNumPositron += 1.0;
                   break;
        // e
        case  11 : fListDataPerThread[indx].fDataGlobal.fNumElectron += 1.0;
                   break;
        // gamma
        case  22 : fListDataPerThread[indx].fDataGlobal.fNumGamma    += 1.0;
                   break;
      }
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
}


// new interface V3 scalar version
void TestEM3ApplicationRP::SteppingActions(GeantTrack &track, GeantTaskData *td) {
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
//    int ilayer = idnode;
    int iabs   = -1;
    if (idvol==fIdAbs) { iabs = 0; }
    else if (idvol==fIdGap) { iabs = 1;}

    const geantphysics::Particle *part = geantphysics::Particle::GetParticleByInternalCode(track.fGVcode);
    int   pdgCode = part->GetPDGCode();
    double charge = part->GetPDGCharge();

    if (iabs>-1) {
      fListDataPerThread[indx].fListDataPerAbsorber[iabs].fEdep     +=  track.fEdep;
      if (charge!=0.0) {
        fListDataPerThread[indx].fListDataPerAbsorber[iabs].fLength +=  track.fStep;
//        std::cout<< track.fBoundary<< " "<<" iabs = "<< iabs<<" "<<std::setprecision(15)<<track.fStep << "  " << track.fTheZPathLenght << "  "<< track.fTheTrueStepLenght << "  "<<track.fPstep<<std::endl;
      }
    }
    if (track.fStatus!=Geant::kNew) { // do not count the creation step
      if (charge==0.0) {
        fListDataPerThread[indx].fDataGlobal.fNumNeutralSteps +=1.0;
      } else {
        fListDataPerThread[indx].fDataGlobal.fNumChargedSteps +=1.0;
      }
    }
    if (track.fStatus==Geant::kNew) { // check secondaries
      switch(pdgCode) {
        // e+
        case -11 : fListDataPerThread[indx].fDataGlobal.fNumPositron += 1.0;
                   break;
        // e
        case  11 : fListDataPerThread[indx].fDataGlobal.fNumElectron += 1.0;
                   break;
        // gamma
        case  22 : fListDataPerThread[indx].fDataGlobal.fNumGamma    += 1.0;
                   break;
      }
    }
//  }

  if (fRunMgr->GetConfig()->fFillTree) {
    MyHit *hit;
    //    int nhits = 0;
//    for (int i = 0; i < npart; i++) {
      // Deposit hits in scintillator
      if (idvol==fIdGap && track.fEdep>0.00002) {
	      hit         = fFactory->NextFree(track.fEvslot, tid);
	      hit->fX     = track.fXpos;
	      hit->fY     = track.fYpos;
	      hit->fZ     = track.fZpos;
	      hit->fEdep  = track.fEdep;
        hit->fTime  = track.fTime;
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


//______________________________________________________________________________
void TestEM3ApplicationRP::FinishRun() {
  double nprim = (double)fRunMgr->GetNprimaries();
//  for (int i=0; i<kNlayers; ++i) {
    for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
      for (int iwthrd=1; iwthrd<fNumWThreads; ++iwthrd) {
        fListDataPerThread[0].fListDataPerAbsorber[iabs].fEdep
        += fListDataPerThread[iwthrd].fListDataPerAbsorber[iabs].fEdep;
        fListDataPerThread[0].fListDataPerAbsorber[iabs].fLength
        +=fListDataPerThread[iwthrd].fListDataPerAbsorber[iabs].fLength;
      }
    }
    for (int iwthrd=1; iwthrd<fNumWThreads; ++iwthrd) {
      fListDataPerThread[0].fDataGlobal.fNumGamma
      += fListDataPerThread[iwthrd].fDataGlobal.fNumGamma;
      fListDataPerThread[0].fDataGlobal.fNumElectron
      += fListDataPerThread[iwthrd].fDataGlobal.fNumElectron;
      fListDataPerThread[0].fDataGlobal.fNumPositron
      += fListDataPerThread[iwthrd].fDataGlobal.fNumPositron;
      fListDataPerThread[0].fDataGlobal.fNumChargedSteps
      += fListDataPerThread[iwthrd].fDataGlobal.fNumChargedSteps;
      fListDataPerThread[0].fDataGlobal.fNumNeutralSteps
      += fListDataPerThread[iwthrd].fDataGlobal.fNumNeutralSteps;
    }
//  }
//  int realNumLayesr = 10;
//  int indxLayes[]   = {2,5,6,7,8,9,10,11,12,13};
  std::cout<< std::endl << std::setw(90) << std::setfill('=') << "" <<"\n";
  std::cout<< "  " << std::setw(30)  << std::setfill('*') << ""
           << "  RESULT OF THE SIMULATION "
           << std::setw(30) <<"" <<"\n";
  std::cout<< std::setw(90) << std::setfill('=') << "" <<std::setfill(' ') <<"\n";
  std::cout<<" [-----"<<std::setw(15) <<" Material "
                      <<std::setw(25) <<" Edep [MeV/#primary] "
                      <<std::setw(40) <<" Charged Track Length [cm] ------]"
                      << std::endl;
//  for (int i=0; i<realNumLayesr; ++i) {
//    int ilayer = indxLayes[i];
//  for (int iabs=0; iabs<kNumAbsorbers; ++iabs) {
    std::cout<< std::left << std::setw(8) << " "
             << std::setw(20) << "Lead:"
             << std::setw(20) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[0].fEdep*1000./nprim
             << std::setw(20) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[0].fLength/nprim
             << std::endl;
    std::cout<< std::left  << std::setw(8) << " "
             << std::setw(20) << "liquidArgon:"
             << std::setw(20) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[1].fEdep*1000./nprim
             << std::setw(20) << std::setprecision(6)
             << fListDataPerThread[0].fListDataPerAbsorber[1].fLength/nprim
             << std::endl;
    std::cout<<" "<<std::setw(90)<<std::setfill('-')<<" "<<std::setfill(' ')<<std::endl;
    std::cout<< std::setw(30) << "     Mean number of Gamma         = " << std::right << std::setw(12)
             << fListDataPerThread[0].fDataGlobal.fNumGamma/nprim << std::endl
             << std::setw(30) << "     Mean number of Electron      = " << std::right << std::setw(12)
             << fListDataPerThread[0].fDataGlobal.fNumElectron/nprim << std::endl
             << std::setw(30) << "     Mean number of Positron      = " << std::right << std::setw(12)
             << fListDataPerThread[0].fDataGlobal.fNumPositron/nprim << std::endl
             << std::setw(30) << "     Mean number of Charged steps = " << std::right << std::setw(12)
             << fListDataPerThread[0].fDataGlobal.fNumChargedSteps/nprim << std::endl
             << std::setw(30) << "     Mean number of Neutral steps = " << std::right << std::setw(12)
             << fListDataPerThread[0].fDataGlobal.fNumNeutralSteps/nprim
             << std::endl;
  std::cout<< std::setw(90) << std::setfill('=') << "" << "\n\n";
}
