
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include <TObject.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TLegend.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStopwatch.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TF1.h>
#include <TObjArray.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TSpectrum.h>
#include <TText.h>
#include <TLatex.h>
#include <TRandom.h>

void plot()
{
  auto c = new TCanvas();
  c->SetGrid();
  TGraphErrors graph1("nudy_U235_ang", "%*lg %lg %lg", "");
  graph1.SetTitle("Angular Distribution;Angle (Deg.);");
  graph1.SetLineWidth(2);
  graph1.DrawClone("E3AL");
  auto c1 = new TCanvas();
  c1->SetGrid();
  TGraphErrors graph2("nudy_U235_ene", "%*lg %lg %lg", "");
  graph2.SetTitle("Energy Distribution;Energy (MeV);");
  graph2.SetLineWidth(2);
  graph2.DrawClone("E3AL");
  auto c2 = new TCanvas();
  c2->SetGrid();
  TGraphErrors graph3("nudy_U235_Sec", "%*lg %lg %lg", "");
  graph3.SetTitle("Fission neutrons;Neutrons;");
  graph3.SetLineWidth(2);
  graph3.DrawClone("E3AL");
}
