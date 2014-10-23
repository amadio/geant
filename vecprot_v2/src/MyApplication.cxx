#include "MyApplication.h"
#include "TGeoNode.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "globals.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TROOT.h"
#include <cassert>

ClassImp(MyApplication)

//______________________________________________________________________________
MyApplication::MyApplication()
              :GeantVApplication(),
               fInitialized(kFALSE),
               fIdGap(0),
               fIdAbs(0),
               fFactory(0)
{
// Ctor..
   GeantFactoryStore *store = GeantFactoryStore::Instance();
   fFactory = store->GetFactory<MyHit>(16); 
   printf("Created factory for MyHit");
   memset(fEdepGap, 0, kNlayers*kMaxThreads*sizeof(Float_t));
   memset(fLengthGap, 0, kNlayers*kMaxThreads*sizeof(Float_t));
   memset(fEdepAbs, 0, kNlayers*kMaxThreads*sizeof(Float_t));
   memset(fLengthAbs, 0, kNlayers*kMaxThreads*sizeof(Float_t));
}

//______________________________________________________________________________
Bool_t MyApplication::Initialize()
{
// Initialize application. Geometry must be loaded.
   if (fInitialized) return kTRUE;
   if (!gGeoManager) {
      Error("Initialize", "Geometry not loaded");
      return kFALSE;
   }
   TGeoVolume *lvGap = gGeoManager->GetVolume("liquidArgon");
   TGeoVolume *lvAbs = gGeoManager->GetVolume("Lead");
   
   if (!lvGap || !lvAbs) {
      Error("Initialize","Logical volumes for gap and absorber not found - do you use the right geometry");
      return kFALSE;
   }
   fIdGap = lvGap->GetNumber();
   fIdAbs = lvAbs->GetNumber();
   fInitialized = kTRUE;
   return kTRUE;
}   
   
//______________________________________________________________________________
void MyApplication::StepManager(Int_t tid, Int_t npart, const GeantTrack_v & tracks)
{
    // Application stepping manager. The thread id has to be used to manage storage
// of hits independently per thread.
   if (!fInitialized) return;     // FOR NOW
   // Loop all tracks, check if they are in the right volume and collect the
   // energy deposit and step length
   TGeoNode const *current;
   Int_t idvol, idnode, ilev;
   for (Int_t i=0; i<npart; i++) {
//      printf("%d=>\n", i);
//      tracks.PrintTrack(i);
      ilev = tracks.fPathV[i]->GetLevel();
      if (ilev<1) continue;
      current = tracks.fPathV[i]->GetCurrentNode();
      if (!current) continue;
      idnode = tracks.fPathV[i]->GetNode(ilev-1)->GetNumber();
      idvol = current->GetVolume()->GetNumber();
      if (idvol == fIdGap) {
//         tracks.PrintTrack(i);
         fEdepGap[idnode-3][tid] += tracks.fEdepV[i];
         fLengthGap[idnode-3][tid] += tracks.fStepV[i];
      } else if (idvol == fIdAbs) {
//         tracks.PrintTrack(i);
         fEdepAbs[idnode-3][tid] += tracks.fEdepV[i];
         fLengthAbs[idnode-3][tid] += tracks.fStepV[i];
      }
   }   
   return;
   MyHit *hit;
   Int_t nhits = 0;
   for (Int_t i=0; i<npart; i++) {
      hit = fFactory->NextFree(tracks.fEvslotV[i]);
      hit->fX = tracks.fXposV[i];
      hit->fY = tracks.fYposV[i];
      hit->fZ = tracks.fZposV[i];
      hit->fEdep = tracks.fEdepV[i];
      hit->fVolId = 0.;
      hit->fDetId = 0.;
//      if (track->path && track->path->GetCurrentNode()) {
//         hit->fVolId = track->path->GetCurrentNode()->GetVolume()->GetNumber();
//         hit->fDetId = track->path->GetCurrentNode()->GetNumber();
//      }
      nhits++;
   }
//   Printf("Thread %d produced %d hits", tid, nhits);
}

//______________________________________________________________________________
void MyApplication::Digitize(Int_t /* event */)
{
// User method to digitize a full event, which is at this stage fully transported
//   printf("======= Statistics for event %d:\n", event);
   Printf("Energy deposit [MeV/primary] and cumulated track length [cm/primary] per layer");
   Printf("================================================================================");
   Double_t nprim = (Double_t)gPropagator->fNprimaries;
   for (Int_t i=0; i<kNlayers; ++i) {
      for (Int_t tid=1; tid<kMaxThreads; ++tid) {
         fEdepGap[i][0] += fEdepGap[i][tid];
         fLengthGap[i][0] += fLengthGap[i][tid];
         fEdepAbs[i][0] += fEdepAbs[i][tid];
         fLengthAbs[i][0] += fLengthAbs[i][tid];
      }
      Printf("Layer %d: Egap=%f   Lgap=%f   Eabs=%f   Labs=%f", i+3, 
             fEdepGap[i][0]*1000./nprim, fLengthGap[i][0]/nprim, fEdepAbs[i][0]*1000./nprim, fLengthAbs[i][0]/nprim);
   }
   Printf("================================================================================");
//   TCanvas *c1 = new TCanvas("Edep", "Energy deposition for ExN03", 700, 800);
   TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("capp");
   if (!c1) return;
   c1->Divide(1,2);
   TVirtualPad *pad = c1->cd(1);
   pad->SetGridx();
   pad->SetGridy();
   pad->SetLogy();
   TH1F *histeg = new TH1F("Edep_gap", "Primary track energy deposition per layer", 12, 0.5, 12.5);
   histeg->SetMarkerColor(kRed);
   histeg->SetMarkerStyle(2);
   histeg->SetStats(kFALSE);
   TH1F *histea = new TH1F("Edep_abs", "Primary track energy deposition per layer in absorber", 12, 0.5, 12.5);
   histea->SetMarkerColor(kBlue);
   histea->SetMarkerStyle(4);
   histea->SetStats(kFALSE);
   for (Int_t i=0; i<10; i++) {
      histeg->SetBinContent(i+3,fEdepGap[i][0]*1000./nprim);
      histea->SetBinContent(i+3,fEdepAbs[i][0]*1000./nprim);
   }
   Double_t minval = TMath::Min(histeg->GetBinContent(histeg->GetMinimumBin()), histea->GetBinContent(histea->GetMinimumBin()));
   minval = TMath::Max(minval, 1.E-5);
   Double_t maxval = TMath::Max(histeg->GetBinContent(histeg->GetMaximumBin()), histea->GetBinContent(histea->GetMaximumBin()));
   histeg->GetXaxis()->SetTitle("Layer");
   histeg->GetYaxis()->SetTitle("Edep per layer [MeV]");
   histeg->GetYaxis()->SetRangeUser(minval-0.1*minval,maxval+0.1*maxval);
   histeg->Draw("P");
   histea->Draw("SAMEP");
//   TCanvas *c2 = new TCanvas("Length", "Length in layers for ExN03", 700, 800);
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
   for (Int_t i=0; i<10; i++) {
      histlg->SetBinContent(i+3,fLengthGap[i][0]/nprim);
      histla->SetBinContent(i+3,fLengthAbs[i][0]/nprim);
   }
   histlg->GetXaxis()->SetTitle("Layer");
   histlg->GetYaxis()->SetTitle("Length per layer");
   minval = TMath::Min(histlg->GetBinContent(histlg->GetMinimumBin()), histla->GetBinContent(histla->GetMinimumBin()));
   minval = TMath::Max(minval, 1.E-5);
   maxval = TMath::Max(histlg->GetBinContent(histlg->GetMaximumBin()), histla->GetBinContent(histla->GetMaximumBin()));
   histlg->GetYaxis()->SetRangeUser(minval-0.1*minval,maxval+0.1*maxval);
   histlg->Draw("P");
   histla->Draw("SAMEP");
}
