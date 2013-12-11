/*#include <TPartIndex.h>
#include <TEXsec.h>
#include <TPXsec.h>
#include <TFile.h>
#include <TH1F.h>*/

void plotxsec(const char* mat, const char* pnam, const char* reac)
{
   gSystem->Load("libXsec");
   TFile *fx = new TFile("xsec.root");
   TPartIndex* tp = fx->Get("PartIndex");
   TEXsec *mate = fx->Get(mat);
   mate->Print();
   //   mate->Dump();
   const Float_t emin = mate->Emin();
   const Float_t emax = mate->Emax();
   const Int_t nbins=100;
   const Float_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
   Float_t *xsec = new Float_t[nbins];
   TH1F *h = new TH1F("h","x-sec",nbins,TMath::Log10(emin),TMath::Log10(emax));
   Double_t en=emin;
   Int_t rcode = tp->ProcIndex(reac);
   Int_t pdg = tp->PDG(pnam);
   Int_t ipart = tp->PartIndex(pnam);
   printf("pdg %d code for %s = %d\n",pdg,reac,rcode);
   for(Int_t i=0; i<nbins; ++i) {
      Float_t xs = mate->XS(ipart,rcode,en);
      printf("xs(%14.7g) = %f\n",en,xs);
      h->SetBinContent(i,xs);
      en*=delta;
   }
   h->Draw();
}
