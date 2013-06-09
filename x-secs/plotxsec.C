/*#include <TPartIndex.h>
#include <TMXsec.h>
#include <TPXsec.h>
#include <TFile.h>
#include <TH1F.h>*/

void plotxsec1(const char* mat, const char* pnam, const char* reac)
{
   gSystem->Load("libXSec");
   TFile *fx = new TFile("xsec.root");
   TPartIndex* tp = fx->Get("PartIndex");
   TMXsec *mate = fx->Get(mat);
   mate->Print();
   //   mate->Dump();
   const Float_t emin = 1e-3;
   const Float_t emax = 1e6;
   const Int_t nbins=100;
   const Float_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
   Float_t *xsec = new Float_t[nbins];
   TH1F *h = new TH1F("h","x-sec",nbins,TMath::Log10(emin),TMath::Log10(emax));
   Double_t en=emin;
   Int_t rcode = tp->Rcode(reac);
   Int_t pdg = tp->PDG(pnam);
   printf("pdg %d code for %s = %d\n",pdg,reac,rcode);
   for(Int_t i=0; i<nbins; ++i) {
      Float_t xs = mate->XS(pdg,rcode,en);
      printf("xs(%14.7g) = %f\n",en,xs);
      h->SetBinContent(i,xs);
      en*=delta;
   }
   h->Draw();
}
