#include "CMSApplication.h"
#include "TGeoNode.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "globals.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <cassert>

ClassImp(CMSApplication)

//______________________________________________________________________________
CMSApplication::CMSApplication()
: GeantVApplication(), fInitialized(kFALSE), fECALMap(), fHCALMap() {
  // Ctor..
  memset(fSensFlags, 0, kNvolumes*sizeof(Bool_t));
  memset(fEdepECAL, 0, kNECALModules * kMaxThreads * sizeof(Float_t));
  memset(fEdepHCAL, 0, kNHCALModules * kMaxThreads * sizeof(Float_t));
  memset(fECALid, 0, kNECALModules*sizeof(Int_t));
  memset(fHCALid, 0, kNHCALModules*sizeof(Int_t));
}

//______________________________________________________________________________
Bool_t CMSApplication::Initialize() {
  // Initialize application. Geometry must be loaded.
  if (fInitialized)
    return kTRUE;
  // Loop unique volume id's
  TGeoVolume *vol;
  TString svol, smat;
  Int_t necal = 0;
  Int_t nhcal = 0;
  for (Int_t ivol=0; ivol<kNvolumes; ++ivol) {
    vol = gGeoManager->GetVolume(ivol);
    svol = vol->GetName();
    // ECAL cells
    if (svol.BeginsWith("EBRY") || svol.BeginsWith("EFRY")) {
      fSensFlags[ivol] = true;
      fECALMap[ivol] = necal;
      fECALid[necal] = ivol;
      necal++;
    }
    // HCAL cells
    if (vol->GetMedium()) smat = vol->GetMaterial()->GetName();
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
void CMSApplication::StepManager(Int_t tid, Int_t npart, const GeantTrack_v &tracks) {
  // Application stepping manager. The thread id has to be used to manage storage
  // of hits independently per thread.
  if (!fInitialized)
    return; // FOR NOW
  // Loop all tracks, check if they are in the right volume and collect the
  // energy deposit and step length
  Int_t ivol;
  Int_t idtype;
  Int_t mod;
  TGeoVolume *vol;
  for (Int_t i = 0; i < npart; i++) {
    vol = tracks.GetVolume(i);
    ivol = vol->GetNumber();
    idtype = 0;
    if (!fSensFlags[ivol]) continue;
    if (vol->GetName()[0] == 'E') idtype = 1;
    else idtype = 2;
    switch (idtype) {
      case 1:
        mod = fECALMap.find(ivol)->second;
        fEdepECAL[mod][tid] += tracks.fEdepV[i];
        break;
      case 2:
        mod = fHCALMap.find(ivol)->second;
        fEdepHCAL[mod][tid] += tracks.fEdepV[i];
        break;
    }
  }
  //   Printf("Thread %d produced %d hits", tid, nhits);
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
    Printf("   volume %s: edep=%f", gGeoManager->GetVolume(fECALid[i])->GetName(), fEdepECAL[i][0] * 1000. / nprim);
  }
  Printf("Energy deposit in HCAL [MeV/primary] ");
  Printf("================================================================================");
  for (Int_t i = 0; i < kNHCALModules; ++i) {
    for (Int_t tid = 1; tid < kMaxThreads; ++tid) {
      fEdepHCAL[i][0] += fEdepHCAL[i][tid];
    }
    Printf("   volume %s: edep=%f", gGeoManager->GetVolume(fHCALid[i])->GetName(), fEdepHCAL[i][0] * 1000. / nprim);
  }
  Printf("================================================================================");
  //   TCanvas *c1 = new TCanvas("Edep", "Energy deposition for CMS", 700, 800);
/*
  TCanvas *c1 = (TCanvas *)gROOT->GetListOfCanvases()->FindObject("capp");
  if (!c1)
    return;
  c1->Divide(1, 2);
  TVirtualPad *pad = c1->cd(1);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetLogy();
  TH1F *histeg = new TH1F("Edep_gap", "Primary track energy deposition per layer", 12, 0.5, 12.5);
  histeg->SetMarkerColor(kRed);
  histeg->SetMarkerStyle(2);
  histeg->SetStats(kFALSE);
  TH1F *histea =
      new TH1F("Edep_abs", "Primary track energy deposition per layer in absorber", 12, 0.5, 12.5);
  histea->SetMarkerColor(kBlue);
  histea->SetMarkerStyle(4);
  histea->SetStats(kFALSE);
  for (Int_t i = 0; i < 10; i++) {
    histeg->SetBinContent(i + 3, fEdepGap[i][0] * 1000. / nprim);
    histea->SetBinContent(i + 3, fEdepAbs[i][0] * 1000. / nprim);
  }
  Double_t minval = TMath::Min(histeg->GetBinContent(histeg->GetMinimumBin()),
                               histea->GetBinContent(histea->GetMinimumBin()));
  minval = TMath::Max(minval, 1.E-5);
  Double_t maxval = TMath::Max(histeg->GetBinContent(histeg->GetMaximumBin()),
                               histea->GetBinContent(histea->GetMaximumBin()));
  histeg->GetXaxis()->SetTitle("Layer");
  histeg->GetYaxis()->SetTitle("Edep per layer [MeV]");
  histeg->GetYaxis()->SetRangeUser(minval - 0.1 * minval, maxval + 0.1 * maxval);
  histeg->Draw("P");
  histea->Draw("SAMEP");
  //   TCanvas *c2 = new TCanvas("Length", "Length in layers for CMS", 700, 800);
  pad = c1->cd(2);
  pad->SetGridx();
  pad->SetGridy();
  pad->SetLogy();
  TH1F *histlg = new TH1F("Len_gap", "Length per layer normalized per primary", 12, 0.5, 12.5);
  histlg->SetMarkerColor(kRed);
  histlg->SetMarkerStyle(2);
  histlg->SetStats(kFALSE);
  TH1F *histla = new TH1F("Len_abs", "Length per layer normalized per primary", 12, 0.5, 12.5);
  histla->SetMarkerColor(kBlue);
  histla->SetMarkerStyle(4);
  histla->SetStats(kFALSE);
  for (Int_t i = 0; i < 10; i++) {
    histlg->SetBinContent(i + 3, fLengthGap[i][0] / nprim);
    histla->SetBinContent(i + 3, fLengthAbs[i][0] / nprim);
  }
  histlg->GetXaxis()->SetTitle("Layer");
  histlg->GetYaxis()->SetTitle("Length per layer");
  minval = TMath::Min(histlg->GetBinContent(histlg->GetMinimumBin()),
                      histla->GetBinContent(histla->GetMinimumBin()));
  minval = TMath::Max(minval, 1.E-5);
  maxval = TMath::Max(histlg->GetBinContent(histlg->GetMaximumBin()),
                      histla->GetBinContent(histla->GetMaximumBin()));
  histlg->GetYaxis()->SetRangeUser(minval - 0.1 * minval, maxval + 0.1 * maxval);
  histlg->Draw("P");
  histla->Draw("SAMEP");
*/
}
