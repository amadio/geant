#ifdef USE_ROOT
#include "LHCbApplication.h"
#ifdef USE_VECGEOM_NAVIGATOR
#include "management/GeoManager.h"
using vecgeom::GeoManager;
#endif
#include "GeantRunManager.h"
#include "GeantEvent.h"
#include "GeantFactoryStore.h"
#include "GeantTrackVec.h"
#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "GeantTaskData.h"
#include "globals.h"
#include "TProfile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include <cassert>

using namespace Geant;

//______________________________________________________________________________
LHCbApplication::LHCbApplication(GeantRunManager *runmgr)
  : GeantVApplication(runmgr), fInitialized(false), fECALMap(), fHCALMap(), fMHist(), fScore(kNoScore), fFluxElec(0),
    fFluxGamma(0), fFluxP(0), fFluxPi(0), fFluxK(0), fEdepElec(0), fEdepGamma(0), fEdepP(0), fEdepPi(0), fEdepK(0), fFactory(0) {
  // Ctor..
  GeantFactoryStore *store = GeantFactoryStore::Instance();
  fFactory = store->GetFactory<MyHit>(16, runmgr->GetNthreadsTotal());
  //
  memset(fSensFlags, 0, kNvolumes * sizeof(bool));
  memset(fEdepECAL, 0, kNECALModules * kMaxThreads * sizeof(float));
  memset(fEdepHCAL, 0, kNHCALModules * kMaxThreads * sizeof(float));
  memset(fECALid, 0, kNECALModules * sizeof(int));
  memset(fHCALid, 0, kNHCALModules * sizeof(int));
  TH1::AddDirectory(false);
  fFluxElec = new TH1F("hFluxElec", "e+/e- flux/primary in ECAL", 50, 0., 2500.);
  fFluxElec->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxElec->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxGamma = new TH1F("hFluxGamma", "Gamma flux/primary in ECAL", 50, 0., 2500.);
  fFluxGamma->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxGamma->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxP = new TH1F("hFluxP", "Proton flux/primary in ECAL", 50, 0., 2500.);
  fFluxP->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxP->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxPi = new TH1F("hFluxPi", "Pion flux/primary in ECAL", 50, 0., 2500.);
  fFluxPi->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxPi->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fFluxK = new TH1F("hFluxK", "Kaon flux/primary in ECAL", 50, 0., 2500.);
  fFluxK->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fFluxK->GetYaxis()->SetTitle("flux [particles/cm^2/primary]");
  fEdepElec = new TH1F("hEdepElec", "Electron energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepElec->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepElec->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepGamma = new TH1F("hEdepGamma", "Gamma energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepGamma->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepGamma->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepP = new TH1F("hEdepP", "Proton energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepP->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepP->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepPi = new TH1F("hEdepPi", "Pion energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepPi->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepPi->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^3/primary]");
  fEdepK = new TH1F("hEdepK", "Kaon energy deposit density/primary in ECAL", 50, 0., 2500.);
  fEdepK->GetXaxis()->SetTitle("Momentum [MeV/c]");
  fEdepK->GetYaxis()->SetTitle("Energy deposit density[MeV/cm^3/primary]");
  TH1::AddDirectory(true);
}

//______________________________________________________________________________
bool LHCbApplication::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return true;
  // Loop unique volume id's
  int nvolumes = fRunMgr->GetNvolumes();
  vector_t<Volume_t const *> &lvolumes = fRunMgr->GetVolumes();
  Printf("Found %d logical volumes", nvolumes);
  const Volume_t *vol;
  TString svol, smat;
  int nvelo = 0;
  int necal = 0;
  int nhcal = 0;
  for (int ivol = 0; ivol < nvolumes; ++ivol) {
    vol = lvolumes[ivol];
    if (!vol) break;
#ifdef USE_VECGEOM_NAVIGATOR
    int idvol = vol->id();
#else
    int idvol = vol->GetNumber();
#endif
    svol = vol->GetName();


    // VELO & TT
    if (svol.BeginsWith("/dd/Geometry/BeforeMagnetRegion/Velo/Sensors/lvVeloRDetectorPU") ||
	svol.BeginsWith("/dd/Geometry/BeforeMagnetRegion/Velo/Sensors/lvVeloPhiDetector") ||
	svol.BeginsWith("/dd/Geometry/BeforeMagnetRegion/Velo/Sensors/lvVeloRDetector") ||
	svol.BeginsWith("/dd/Geometry/BeforeMagnetRegion/TT/Modules/lvSensor") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/IT/Ladder/lvLong") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/IT/Ladder/lvShort") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/OT/Modules/lvS3ModuleA") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/OT/Modules/lvS1ModuleA") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/OT/Modules/lvLModuleA") ||
	svol.BeginsWith("/dd/Geometry/AfterMagnetRegion/T/OT/Modules/lvS2ModuleA"))
      {
      fSensFlags[idvol] = true;
      fVELOMap[idvol] = nvelo;
      fVELOid[necal] = idvol;
      nvelo++;
    }
    
    // ECAL cells
    if (svol.BeginsWith("/dd/Geometry/DownstreamRegion/Ecal/Modules/InnCell") ||
	svol.BeginsWith("/dd/Geometry/DownstreamRegion/Ecal/Modules/MidCell") ||
	svol.BeginsWith("/dd/Geometry/DownstreamRegion/Ecal/Modules/OutCell")) {
      fSensFlags[idvol] = true;
      fECALMap[idvol] = necal;
      fECALid[necal] = idvol;
      necal++;
    }
    
    // HCAL cells
    if (svol.BeginsWith("/dd/Geometry/DownstreamRegion/Hcal/Cells/lvHcalInnCellUpScTile") ||
	svol.BeginsWith("/dd/Geometry/DownstreamRegion/Hcal/Cells/lvHcalInnCellLowScTile") || 
	svol.BeginsWith("/dd/Geometry/DownstreamRegion/Hcal/Cells/lvHcalOutCellScTile")) {
      fSensFlags[idvol] = true;
      fHCALMap[idvol] = nhcal;
      fHCALid[nhcal] = idvol;
      nhcal++;
    }
  }
  
  Printf("=== LHCbApplication::Initialize: necal=%d  nhcal=%d", necal, nhcal);
  fInitialized = true;  
  return true;
}

//______________________________________________________________________________
void LHCbApplication::StepManager(int npart, const GeantTrack_v &tracks, GeantTaskData *td) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  GeantPropagator *propagator = td->fPropagator;
  int tid = td->fTid;
  if ((!fInitialized) || (fScore == kNoScore))
    return;
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  int ivol;
  int idtype;
  int mod;
  Volume_t const *vol;

  if (propagator->fConfig->fFillTree) {
    for (int itr = 0; itr < npart; itr++) {
      vol = tracks.GetVolume(itr);
#ifdef USE_VECGEOM_NAVIGATOR
      ivol = vol->id();
#else
      ivol = vol->GetNumber();
#endif

      idtype = 0;
      if (fSensFlags[ivol]) {

	if (vol->GetName()[30] == 'E')
	  idtype = 1;
	else if (vol->GetName()[30] == 'H')
	  idtype = 2;
      
	switch (idtype) {
	case 1:
	  mod = fECALMap.find(ivol)->second;
	  fEdepECAL[mod][tid] += tracks.fEdepV[itr];
	  break;
	case 2:
	  mod = fHCALMap.find(ivol)->second;
	  fEdepHCAL[mod][tid] += tracks.fEdepV[itr];
	  break;
	}
      
    
	MyHit *hit;
	// Deposit hits
	//      if (idtype==1) {
	//	Printf("hit at z %f id %i", tracks.fZposV[itr], idtype);
        hit = fFactory->NextFree(tracks.fEvslotV[itr], tid);
        hit->fX = tracks.fXposV[itr];
        hit->fY = tracks.fYposV[itr];
        hit->fZ = tracks.fZposV[itr];
        hit->fEdep = 1000*tracks.fEdepV[itr];
        hit->fTime = tracks.fTimeV[itr];
        hit->fEvent = tracks.fEventV[itr];
        hit->fTrack = tracks.fParticleV[itr];
        hit->fVolId = ivol;
        hit->fDetId = idtype;
	//}

	// Score in ECAL
	if (idtype == 1) {
	  // Add scored entity
      
	  if (propagator->fNthreads > 1)
	    fMHist.lock();
      
	  double capacity = 1.;
    
	  if (fabs(tracks.fPDGV[itr]) == 11) {
	    fFluxElec->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
	    fEdepElec->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
	  } else if (tracks.fPDGV[itr] == 22 || tracks.fPDGV[itr] == 0) {
	    //        std::cout << tracks.GetVolume(itr)->GetName() << " " << tracks.fStepV[itr] << std::endl;
	    fFluxGamma->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
	    fEdepGamma->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
	  } else if (tracks.fPDGV[itr] == 2212) {
	    fFluxP->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
	    fEdepP->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
	  } else if (fabs(tracks.fPDGV[itr]) == 211) {
	    fFluxPi->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
	    fEdepPi->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
	  } else if (fabs(tracks.fPDGV[itr]) == 321) {
	    fFluxK->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
	    fEdepK->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
	  }
	  if (propagator->fNthreads > 1)
	    fMHist.unlock();      
	}
      }
    }
  }
}

//______________________________________________________________________________
void LHCbApplication::FinishEvent(int evt, int slot) {
  // User method to digitize a full event, which is at this stage fully transported
  //   Printf("======= Statistics for event %d:\n", event);
  GeantEvent *event = fRunMgr->GetEvent(slot);
  Printf("Energy deposit for event %d in ECAL [MeV/primary] ", evt);
  Printf("================================================================================");
  double nprim = (double)(event->GetNtracks());
  for (int i = 0; i < kNECALModules; ++i) {
    for (int tid = 1; tid < kMaxThreads; ++tid) {
      fEdepECAL[i][0] += fEdepECAL[i][tid];
    }
#ifdef USE_VECGEOM_NAVIGATOR
    Printf("   volume %s: edep=%f", GeoManager::Instance().FindLogicalVolume(fECALid[i])->GetName(),
           fEdepECAL[i][0] * 1000. / nprim);
#else
    Printf("   volume %s: edep=%f", gGeoManager->GetVolume(fECALid[i])->GetName(), fEdepECAL[i][0] * 1000. / nprim);
#endif
  }
  Printf("Energy deposit in HCAL [MeV/primary] ");
  Printf("================================================================================");
  for (int i = 0; i < kNHCALModules; ++i) {
    for (int tid = 1; tid < kMaxThreads; ++tid) {
      fEdepHCAL[i][0] += fEdepHCAL[i][tid];
    }
#ifdef USE_VECGEOM_NAVIGATOR
    Printf("   volume %s: edep=%f", GeoManager::Instance().FindLogicalVolume(fHCALid[i])->GetName(),
           fEdepHCAL[i][0] * 1000. / nprim);
#else
    Printf("   volume %s: edep=%f", gGeoManager->GetVolume(fHCALid[i])->GetName(), fEdepHCAL[i][0] * 1000. / nprim);
#endif
  }
  Printf("================================================================================");
}

//______________________________________________________________________________
void LHCbApplication::FinishRun() {
  if (fScore == kNoScore)
    return;
  TCanvas *c1 = new TCanvas("LHCb test flux", "Simple scoring in LHCb geometry", 700, 1200);
  double norm = 1. / fRunMgr->GetNprimaries();
  TVirtualPad *pad;
  TFile *f = TFile::Open("ScoreECAL.root", "RECREATE");
  c1->Divide(2, 3);
  pad = c1->cd(1);
  pad->SetLogx();
  pad->SetLogy();
  fFluxElec->Sumw2();
  fFluxElec->Scale(norm);
  fFluxElec->Draw("9");
  pad = c1->cd(2);
  pad->SetLogx();
  pad->SetLogy();
  fFluxGamma->Sumw2();
  fFluxGamma->Scale(norm);
  fFluxGamma->Draw("9");
  pad = c1->cd(3);
  pad->SetLogx();
  pad->SetLogy();
  fFluxP->Sumw2();
  fFluxP->Scale(norm);
  fFluxP->Draw("9");
  pad = c1->cd(4);
  pad->SetLogx();
  pad->SetLogy();
  fFluxPi->Sumw2();
  fFluxPi->Scale(norm);
  fFluxPi->Draw("9");
  pad = c1->cd(5);
  pad->SetLogx();
  pad->SetLogy();
  fFluxK->Sumw2();
  fFluxK->Scale(norm);
  fFluxK->Draw("9");
  fFluxElec->Write();
  fFluxGamma->Write();
  fFluxP->Write();
  fFluxPi->Write();
  fFluxK->Write();

  TCanvas *c2 = new TCanvas("LHCb test edep", "Simple scoring in LHCb geometry", 700, 1200);
  c2->Divide(2, 3);
  pad = c2->cd(1);
  pad->SetLogx();
  pad->SetLogy();
  fEdepElec->Sumw2();
  fEdepElec->Scale(norm);
  fEdepElec->Draw("9");
  pad = c2->cd(2);
  pad->SetLogx();
  pad->SetLogy();
  fEdepP->Sumw2();
  fEdepP->Scale(norm);
  fEdepP->Draw("9");
  pad = c2->cd(3);
  pad->SetLogx();
  pad->SetLogy();
  fEdepPi->Sumw2();
  fEdepPi->Scale(norm);
  fEdepPi->Draw("9");
  pad = c2->cd(4);
  pad->SetLogx();
  pad->SetLogy();
  fEdepK->Sumw2();
  fEdepK->Scale(norm);
  fEdepK->Draw("9");
  pad = c2->cd(5);
  pad->SetLogx();
  pad->SetLogy();
  fEdepGamma->Sumw2();
  fEdepGamma->Scale(norm);
  fEdepGamma->Draw("9");
  fEdepElec->Write();
  fEdepGamma->Write();
  fEdepP->Write();
  fEdepPi->Write();
  fEdepK->Write();

  // Close file
  f->Close();
}
#endif
