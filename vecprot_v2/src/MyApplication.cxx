#include "MyApplication.h"
#include "TGeoBranchArray.h"
#include "TGeoNode.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "globals.h"
#include "TH1.h"
#include "TCanvas.h"

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
   memset(fEdepGap, 0, kNlayers*sizeof(Float_t));
   memset(fLengthGap, 0, kNlayers*sizeof(Float_t));
   memset(fEdepAbs, 0, kNlayers*sizeof(Float_t));
   memset(fLengthAbs, 0, kNlayers*sizeof(Float_t));
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
   TGeoNode *current;
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
         fEdepGap[idnode-3] += tracks.fEdepV[i];
         fLengthGap[idnode-3] += tracks.fStepV[i];
      } else if (idvol == fIdAbs) {
//         tracks.PrintTrack(i);
         fEdepAbs[idnode-3] += tracks.fEdepV[i];
         fLengthAbs[idnode-3] += tracks.fStepV[i];
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
void MyApplication::Digitize(Int_t event)
{
// User method to digitize a full event, which is at this stage fully transported
   printf("======= Statistics for event %d:\n", event);
   for (Int_t i=0; i<kNlayers; ++i) {
      printf("Layer %d: Egap=%f   Lgap=%f   Eabs=%f   Labs=%f\n", i+3, 
             fEdepGap[i], fLengthGap[i], fEdepAbs[i], fLengthAbs[i]);
   }
//   memset(fEdepGap, 0, kNlayers*sizeof(Float_t));
//   memset(fLengthGap, 0, kNlayers*sizeof(Float_t));
//   memset(fEdepAbs, 0, kNlayers*sizeof(Float_t));
//   memset(fLengthAbs, 0, kNlayers*sizeof(Float_t));
   TCanvas *c1 = new TCanvas("Edep", "Energy deposition for ExN03", 700, 800);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLogy();
   TH1F *histeg = new TH1F("Edep_gap", "Primary track energy deposition per layer", 10, -0.5, 12.5);
   histeg->SetMarkerColor(kRed);
   histeg->SetMarkerStyle(2);
   histeg->SetStats(kFALSE);
   TH1F *histea = new TH1F("Edep_abs", "Primary track energy deposition per layer in absorber", 10, -0.5, 12.5);
   histea->SetMarkerColor(kBlue);
   histea->SetMarkerStyle(4);
   histea->SetStats(kFALSE);
   for (Int_t i=0; i<10; i++) {
      histeg->SetBinContent(i+3,fEdepGap[i]*1000./(Double_t)gPropagator->fNprimaries);
      histea->SetBinContent(i+3,fEdepAbs[i]*1000./(Double_t)gPropagator->fNprimaries);
   }
   Double_t minval = TMath::Min(histeg->GetBinContent(histeg->GetMinimumBin()), histea->GetBinContent(histea->GetMinimumBin()));
   minval = TMath::Max(minval, 1.E-5);
   Double_t maxval = TMath::Max(histeg->GetBinContent(histeg->GetMaximumBin()), histea->GetBinContent(histea->GetMaximumBin()));
   histeg->GetXaxis()->SetTitle("Layer");
   histeg->GetYaxis()->SetTitle("Edep per layer [MeV]");
   histeg->GetYaxis()->SetRangeUser(minval-0.1*minval,maxval+0.1*maxval);
   histeg->Draw("P");
   histea->Draw("SAMEP");
   TCanvas *c2 = new TCanvas("Length", "Length in layers for ExN03", 700, 800);
   c2->SetGridx();
   c2->SetGridy();
   c2->SetLogy();
   TH1F *histlg = new TH1F("Len_gap", "Length per layer normalized per primary", 10, -0.5, 12.5);
   histlg->SetMarkerColor(kRed);
   histlg->SetMarkerStyle(2);
   histlg->SetStats(kFALSE);
   TH1F *histla = new TH1F("Len_abs", "Length per layer normalized per primary", 10, -0.5, 12.5);
   histla->SetMarkerColor(kBlue);
   histla->SetMarkerStyle(4);
   histla->SetStats(kFALSE);
   for (Int_t i=0; i<10; i++) {
      histlg->SetBinContent(i+3,fLengthGap[i]/(Double_t)gPropagator->fNprimaries);
      histla->SetBinContent(i+3,fLengthAbs[i]/(Double_t)gPropagator->fNprimaries);
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
