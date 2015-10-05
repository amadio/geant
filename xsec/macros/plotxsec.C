#include "TGeoElement.h"
#include "TGeoMaterial.h"
#include "TPartIndex.h"
#include "TSystem.h"
#include "TEXsec.h"
#include "TPXsec.h"
#include "TFile.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TCanvas.h"

#include "string"

void plotxsec(const char* mat="Fe", const char* pnam="proton", const char* reac="inElastic")
{
   const char *fname = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
   gSystem->Load("libXsec");
   TFile *fx = new TFile(fname,"read");
   fx->Get("PartIndex");
   TPartIndex *tp = TPartIndex::I();
   TEXsec *mate = nullptr;
   fx->GetObject(mat,mate);
   //   mate->Print();
   //   mate->Dump();
   const Float_t emin = mate->Emin();
   const Float_t emax = mate->Emax();
   const Int_t nbins=100;
   double *ven = new double[nbins];
   double *vxs = new double[nbins];
   const Float_t delta = exp(log(emax/emin)/(nbins-1));
   Float_t *xsec = new Float_t[nbins];
   std::string title(reac);
   title += " cross section of ";
   title += pnam;
   title += " on ";
   title += mat;
   Double_t en=emin;
   Int_t rcode = tp->ProcIndex(reac);
   Int_t pdg = tp->PDG(pnam);
   Int_t ipart = tp->PartIndex(pnam);
   printf("pdg %d code for %s = %d\n",pdg,reac,rcode);
   for(Int_t i=0; i<nbins; ++i) {
      Float_t xs = mate->XS(ipart,rcode,en);
      printf("xs(%14.7g) = %f\n",en,xs);
      ven[i] = en;
      vxs[i] = xs;
      en*=delta;
   }

   TGraph *g = new TGraph(nbins, ven, vxs);

   TCanvas *c1 = new TCanvas();
   c1->SetLogx();
   g->SetTitle(title.c_str());
   g->GetXaxis()->SetTitle("E [GeV]");
   g->GetYaxis()->SetTitle("#sigma [barn]");
   g->SetMarkerStyle(20);
   g->Draw("apc");
}
