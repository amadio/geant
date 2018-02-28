#include "TError.h"
#include "TFile.h"
#include "TLegend.h"

void DrawComparison(TH1F *hg4, TH1F *hgv, TVirtualPad *pad);

void compareG4_GV()
{
  // Draw comparison plot for the particle flux end energy deposit in the CMS
  // ECAL, for Geant4 anf GeantV
  TFile *fg4 = TFile::Open("hists.root");
  if (!fg4) {
    ::Error("compareG4_GV", "File hists.root not found");
    return;
  }
  TFile *fgv = TFile::Open("ScoreECAL.root");
  if (!fgv) {
    ::Error("compareG4_GV", "File ScoreECAL.root not found");
    return;
  }
  TH1F *hG4FluxElec = (TH1F *)fg4->Get("hFluxElec");
  hG4FluxElec->SetMarkerColor(1);
  hG4FluxElec->SetMarkerStyle(4);
  TH1F *hG4FluxGamma = (TH1F *)fg4->Get("hFluxGamma");
  hG4FluxGamma->SetMarkerColor(1);
  hG4FluxGamma->SetMarkerStyle(4);
  TH1F *hG4FluxP = (TH1F *)fg4->Get("hFluxP");
  hG4FluxP->SetMarkerColor(1);
  hG4FluxP->SetMarkerStyle(4);
  TH1F *hG4FluxPi = (TH1F *)fg4->Get("hFluxPi");
  hG4FluxPi->SetMarkerColor(1);
  hG4FluxPi->SetMarkerStyle(4);
  TH1F *hG4FluxK = (TH1F *)fg4->Get("hFluxK");
  hG4FluxK->SetMarkerColor(1);
  hG4FluxK->SetMarkerStyle(4);
  TH1F *hGVFluxElec = (TH1F *)fgv->Get("hFluxElec");
  hGVFluxElec->SetMarkerColor(kRed);
  hGVFluxElec->SetMarkerStyle(20);
  TH1F *hGVFluxGamma = (TH1F *)fgv->Get("hFluxGamma");
  hGVFluxGamma->SetMarkerColor(kRed);
  hGVFluxGamma->SetMarkerStyle(20);
  TH1F *hGVFluxP = (TH1F *)fgv->Get("hFluxP");
  hGVFluxP->SetMarkerColor(kRed);
  hGVFluxP->SetMarkerStyle(20);
  TH1F *hGVFluxPi = (TH1F *)fgv->Get("hFluxPi");
  hGVFluxPi->SetMarkerColor(kRed);
  hGVFluxPi->SetMarkerStyle(20);
  TH1F *hGVFluxK = (TH1F *)fgv->Get("hFluxK");
  hGVFluxK->SetMarkerColor(kRed);
  hGVFluxK->SetMarkerStyle(20);

  TVirtualPad *pad;
  TCanvas *c1 = new TCanvas("CMS test", "Flux scoring comparison in CMS ECAL", 1200, 800);
  c1->Divide(2, 3);
  pad = c1->cd(1);
  DrawComparison(hG4FluxElec, hGVFluxElec, pad);

  pad = c1->cd(2);
  DrawComparison(hG4FluxGamma, hGVFluxGamma, pad);

  pad = c1->cd(3);
  DrawComparison(hG4FluxP, hGVFluxP, pad);

  pad = c1->cd(4);
  DrawComparison(hG4FluxPi, hGVFluxPi, pad);

  pad = c1->cd(5);
  DrawComparison(hG4FluxK, hGVFluxK, pad);

  // Edep
  TH1F *hG4EdepElec = (TH1F *)fg4->Get("hEdepElec");
  hG4EdepElec->SetMarkerColor(1);
  hG4EdepElec->SetMarkerStyle(4);
  TH1F *hG4EdepGamma = (TH1F *)fg4->Get("hEdepGamma");
  hG4EdepGamma->SetMarkerColor(1);
  hG4EdepGamma->SetMarkerStyle(4);
  TH1F *hG4EdepP = (TH1F *)fg4->Get("hEdepP");
  hG4EdepP->SetMarkerColor(1);
  hG4EdepP->SetMarkerStyle(4);
  TH1F *hG4EdepPi = (TH1F *)fg4->Get("hEdepPi");
  hG4EdepPi->SetMarkerColor(1);
  hG4EdepPi->SetMarkerStyle(4);
  TH1F *hG4EdepK = (TH1F *)fg4->Get("hEdepK");
  hG4EdepK->SetMarkerColor(1);
  hG4EdepK->SetMarkerStyle(4);
  TH1F *hGVEdepElec = (TH1F *)fgv->Get("hEdepElec");
  hGVEdepElec->SetMarkerColor(kRed);
  hGVEdepElec->SetMarkerStyle(20);
  TH1F *hGVEdepGamma = (TH1F *)fgv->Get("hEdepGamma");
  hGVEdepGamma->SetMarkerColor(kRed);
  hGVEdepGamma->SetMarkerStyle(20);
  TH1F *hGVEdepP = (TH1F *)fgv->Get("hEdepP");
  hGVEdepP->SetMarkerColor(kRed);
  hGVEdepP->SetMarkerStyle(20);
  TH1F *hGVEdepPi = (TH1F *)fgv->Get("hEdepPi");
  hGVEdepPi->SetMarkerColor(kRed);
  hGVEdepPi->SetMarkerStyle(20);
  TH1F *hGVEdepK = (TH1F *)fgv->Get("hEdepK");
  hGVEdepK->SetMarkerColor(kRed);
  hGVEdepK->SetMarkerStyle(20);

  TCanvas *c2 = new TCanvas("CMS test edep", "Edep scoring comparison in CMS ECAL", 1600, 1200);
  c2->Divide(2, 2);
  pad = c2->cd(1);
  DrawComparison(hG4EdepElec, hGVEdepElec, pad);

  //  pad = c2->cd(2);
  //  DrawComparison(hG4EdepGamma, hGVEdepGamma, pad);

  pad = c2->cd(2);
  DrawComparison(hG4EdepP, hGVEdepP, pad);

  pad = c2->cd(3);
  DrawComparison(hG4EdepPi, hGVEdepPi, pad);

  pad = c2->cd(4);
  DrawComparison(hG4EdepK, hGVEdepK, pad);

  TCanvas *c3 = new TCanvas("CMS energy deposit protons", "Energy deposit density comparison in CMS ECAL", 1280, 1024);
  hG4EdepP->GetYaxis()->SetTitle("Energy deposit density [MeV/cm^{3}/primary]");
  DrawComparison(hG4EdepP, hGVEdepP, c3->cd());

  // Speed comparison plot
  Double_t x[5] = {5, 10, 20, 50, 100};
  //  Double_t y[5] = {1.23, 1.34, 1.34, 1.34, 1.36};
  Double_t y[5] = {1.99, 1.99, 1.81, 1.89, 2.36};
  TGraph *gr    = new TGraph(5, x, y);
  gr->SetTitle("Single thread performance comparison");
  gr->GetXaxis()->SetTitle("# events");
  gr->GetYaxis()->SetTitle("Time_{Geant4}/Time_{GeantV}");
  gr->GetYaxis()->SetRangeUser(1.5, 2.5);
  gr->GetYaxis()->SetNdivisions(505);
  TCanvas *c4 =
      new TCanvas("Simulation time Geant4/GeantV", "Energy deposit density comparison in CMS ECAL", 1600, 1200);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerSize(1.4);
  gr->SetLineColor(kRed);
  gr->SetLineWidth(2);
  c4->SetGridy();
  gr->Draw("AL*");
  gr->SetMarkerStyle(20);

  c1->SaveAs("fluxECAL.pdf");
  c2->SaveAs("edepECAL.pdf");
  c3->SaveAs("edepProtonECAL.pdf");
  c4->SaveAs("perfG4_GV.pdf");
}

void DrawComparison(TH1F *hg4, TH1F *hgv, TVirtualPad *pad)
{
  TLegend *legend = new TLegend(0.69, 0.75, 0.89, 0.88);
  legend->SetLineColor(0);
  legend->AddEntry(hg4, "Geant4", "lpe");
  legend->AddEntry(hgv, "GeantV", "lpe");
  pad->cd();
  TPad *padtop = new TPad("top", "", 0, 0.3, 1, 1);
  padtop->SetBottomMargin(0);
  padtop->SetLogx();
  padtop->SetLogy();
  //  padtop->SetGridx();
  //  padtop->SetGridy();
  padtop->Draw();
  padtop->cd();
  hg4->SetStats(0);
  hg4->GetYaxis()->SetTitleSize(0.04);
  hg4->GetYaxis()->SetTitleOffset(1.1);
  hg4->GetYaxis()->SetLabelSize(0.04);
  hg4->GetYaxis()->SetLabelOffset(0.005);
  hg4->GetXaxis()->SetRangeUser(50., 2500.);
  hgv->GetXaxis()->SetRangeUser(50., 2500.);
  hg4->Draw("9");
  hgv->Draw("9SAME");
  legend->Draw();

  pad->cd();
  TPad *padbottom = new TPad("bottom", "", 0, 0.05, 1, 0.3);
  padbottom->SetTopMargin(0);
  padbottom->SetBottomMargin(0.25);
  //  padbottom->SetGridx();
  padbottom->SetGridy();
  padbottom->SetLogx();
  padbottom->Draw();
  padbottom->cd();
  TString hname = hgv->GetName();
  hname += "Rap";
  TH1F *hRap = (TH1F *)hgv->Clone(hname);
  hRap->SetTitle("");
  hRap->GetYaxis()->SetTitle("ratio ");
  hRap->GetYaxis()->SetNdivisions(505);
  hRap->GetYaxis()->SetTitleSize(0.12);
  hRap->GetYaxis()->SetTitleOffset(0.37);
  hRap->GetYaxis()->SetLabelSize(0.11);
  hRap->GetYaxis()->SetLabelOffset(0.005);
  hRap->GetXaxis()->SetTitleSize(0.12);
  hRap->GetXaxis()->SetTitleOffset(0.8);
  hRap->GetXaxis()->SetLabelSize(0.12);
  hRap->GetXaxis()->SetLabelOffset(0.002);
  hRap->GetXaxis()->CenterTitle();
  hRap->Sumw2();
  hRap->SetMarkerColor(kBlack);
  hRap->SetMarkerSize(0.8);
  hRap->SetMarkerStyle(21);
  hRap->SetStats(0);
  hRap->Divide(hg4);
  hRap->Draw("ep");
}
