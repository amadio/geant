#include "CMSApplication.h"
#include "ExN03Application.h"
#ifdef USE_VECGEOM_NAVIGATOR
 #include "management/GeoManager.h"
 using vecgeom::GeoManager;
#endif
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantTaskData.h"
#include "globals.h"
#include "TProfile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include <cassert>

ClassImp(CMSApplication)

//______________________________________________________________________________
CMSApplication::CMSApplication()
: GeantVApplication(), fInitialized(kFALSE), fECALMap(), fHCALMap(), fMHist(), fScore(kNoScore), fFluxElec(0),
  fFluxGamma(0), fFluxP(0), fFluxPi(0), fFluxK(0), fEdepElec(0), fEdepGamma(0), fEdepP(0), fEdepPi(0), fEdepK(0) {
  // Ctor..
  memset(fSensFlags, 0, kNvolumes * sizeof(Bool_t));
  memset(fEdepECAL, 0, kNECALModules * kMaxThreads * sizeof(Float_t));
  memset(fEdepHCAL, 0, kNHCALModules * kMaxThreads * sizeof(Float_t));
  memset(fECALid, 0, kNECALModules * sizeof(Int_t));
  memset(fHCALid, 0, kNHCALModules * sizeof(Int_t));
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
Bool_t CMSApplication::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return kTRUE;
  // Loop unique volume id's
  Volume_t *vol;
  TString svol, smat;
  Int_t necal = 0;
  Int_t nhcal = 0;
  for (Int_t ivol = 0; ivol < kNvolumes; ++ivol) {
#ifdef USE_VECGEOM_NAVIGATOR
    vol = GeoManager::Instance().FindLogicalVolume(ivol);
#else
    vol = gGeoManager->GetVolume(ivol);
#endif
    svol = vol->GetName();
    // ECAL cells
    if (svol.BeginsWith("EBRY") || svol.BeginsWith("EFRY")) {
      fSensFlags[ivol] = true;
      fECALMap[ivol] = necal;
      fECALid[necal] = ivol;
      necal++;
    }
// HCAL cells
#ifdef USE_VECGEOM_NAVIGATOR
    //    cout << __func__ << "::vol " << vol->GetName() << endl;
    if (vol->getTrackingMediumPtr())
      smat = ((vecgeom::Medium *)vol->getTrackingMediumPtr())->GetMaterial()->GetName();
#else
    if (vol->GetMedium())
      smat = vol->GetMaterial()->GetName();
#endif
    if (smat == "Scintillator") {
      fSensFlags[ivol] = true;
      fHCALMap[ivol] = nhcal;
      fHCALid[nhcal] = ivol;
      nhcal++;
    }
  }
  Printf("=== CMSApplication::Initialize: necal=%d  nhcal=%d", necal, nhcal);
  fInitialized = kTRUE;
  return kTRUE;
}

//______________________________________________________________________________
void CMSApplication::StepManager(Int_t npart, const GeantTrack_v &tracks, GeantTaskData *td) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  static GeantPropagator *propagator = GeantPropagator::Instance();
  Int_t tid = td->fTid;
  if ((!fInitialized) || (fScore == kNoScore))
    return;
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  Int_t ivol;
  Int_t idtype;
  Int_t mod;
  Volume_t *vol;
  for (Int_t itr = 0; itr < npart; itr++) {
    vol = tracks.GetVolume(itr);
#ifdef USE_VECGEOM_NAVIGATOR
    ivol = vol->id();
#else
    ivol = vol->GetNumber();
#endif
    idtype = 0;
    if (fSensFlags[ivol]) {
      if (vol->GetName()[0] == 'E')
        idtype = 1;
      else
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
    }
    // Score in ECAL
    if (idtype == 1) {
      // Add scored entity
      if (propagator->fNthreads > 1)
        fMHist.lock();
      Double_t capacity = 0.;
#ifdef USE_VECGEOM_NAVIGATOR
      capacity = 1.;
#else
      capacity = vol->GetShape()->Capacity();
#endif
      if (fabs(tracks.fPDGV[itr]) == 11) {
        fFluxElec->Fill(1000. * tracks.fPV[itr], tracks.fStepV[itr] / capacity);
        fEdepElec->Fill(1000. * tracks.fPV[itr], 1000. * tracks.fEdepV[itr] / capacity);
      } else if (tracks.fPDGV[itr] == 22) {
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

//______________________________________________________________________________
void CMSApplication::Digitize(Int_t /* event */) {
  // User method to digitize a full event, which is at this stage fully transported
  //   printf("======= Statistics for event %d:\n", event);
  Printf("Energy deposit in ECAL [MeV/primary] ");
  Printf("================================================================================");
  Double_t nprim = (Double_t)gPropagator->fNprimaries;
  for (Int_t i = 0; i < kNECALModules; ++i) {
    for (Int_t tid = 1; tid < kMaxThreads; ++tid) {
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
  for (Int_t i = 0; i < kNHCALModules; ++i) {
    for (Int_t tid = 1; tid < kMaxThreads; ++tid) {
      fEdepHCAL[i][0] += fEdepHCAL[i][tid];
    }
#ifdef USE_VECGEOM_NAVIGATOR
    Printf("   volume %s: edep=%f", GeoManager::Instance().FindLogicalVolume(fECALid[i])->GetName(),
           fEdepHCAL[i][0] * 1000. / nprim);
#else
    Printf("   volume %s: edep=%f", gGeoManager->GetVolume(fHCALid[i])->GetName(), fEdepHCAL[i][0] * 1000. / nprim);
#endif
  }
  Printf("================================================================================");
}

//______________________________________________________________________________
void CMSApplication::FinishRun() {
  if (fScore == kNoScore)
    return;
  TCanvas *c1 = new TCanvas("CMS test flux", "Simple scoring in CMS geometry", 700, 1200);
  Double_t norm = 1. / GeantPropagator::Instance()->fNprimaries.load();
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

  TCanvas *c2 = new TCanvas("CMS test edep", "Simple scoring in CMS geometry", 700, 1200);
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
