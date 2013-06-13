#include <TMath.h>
#include <TPartIndex.h>
#include <TMXsec.h>
#include <TH1F.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TMultiGraph.h>

Int_t rsampling(const char *el="O",const char *part="proton",Int_t nrep=100000) 
{
   const EColor col[14] = {kBlack,kGray,kRed,kGreen,kBlue, kMagenta, kCyan,
			   kOrange, kSpring, kTeal, kAzure, kViolet, kPink };

   gSystem->Load("libXsec");
   Double_t emin = 1e-3;
   Double_t emax = 1e6;
   Int_t nbins = 100;
   Double_t edelta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
   // Sampling of reactions
   TFile *f = new TFile("xsec.root");
   TPartIndex *tp = (TPartIndex*) f->Get("PartIndex");
   Double_t *e = new Double_t[nbins];
   Double_t *x = new Double_t[nbins];
   TH1F **histos = new TH1F*[TPartIndex::I()->NReac()];

   TMultiGraph *tmg= new TMultiGraph("G5","G5");
   TGraph *gr=0;

   TMXsec * o = (TMXsec*)f->Get(el);
   Double_t en=emin;
   for(Int_t nb=0; nb<nbins; ++nb) {
      e[nb] = TMath::Log10(en);
      en*=edelta;
   }
   Int_t pdg = TPartIndex::I()->PDG(part);
   for(Int_t i=0; i<tp->NReac()-1; ++i) {
      Int_t rcode = TPartIndex::I()->ProcCode(i);
      if(o->XSPDG(pdg,rcode,1)>=0) {
	 en=emin;
	 for(Int_t nb=0; nb<nbins; ++nb) {
	    x[nb] = o->XSPDG(pdg, rcode, en);
	    en*=edelta;
	 }
	 gr = new TGraph(nbins,e,x);
	 gr->SetTitle(TPartIndex::I()->ProcNameCode(rcode));
	 gr->SetTitle(TPartIndex::I()->ProcNameCode(rcode));
	 gr->SetLineColor(col[i%14]);
	 tmg->Add(gr);
      }
   }
   
   for(Int_t i=0; i<tp->NReac(); ++i) {
      printf("%d",i);
      histos[i]= new TH1F(tp->ProcNameIndex(i),tp->ProcNameIndex(i),nbins,
			  TMath::Log10(emin),TMath::Log10(emax));
   }
   Int_t ibin=nbins;
   en=emin;
   while(ibin--) {
      Double_t hnorm = o->XSPDG(pdg,999,en)/nrep;
      Double_t len = TMath::Log10(en);
      Int_t irep = nrep;
      while(irep--) {
	 Int_t reac = o->SampleReacPDG(pdg, en);
	 histos[reac]->Fill(len,hnorm);
      }
      en*=edelta;
   }
   Double_t ymax = -1;
   for(Int_t i=0; i<tp->NReac(); ++i)
      ymax = histos[i]->GetMaximum()>ymax?histos[i]->GetMaximum():ymax;
   
   histos[0]->SetMaximum(1.5*ymax);
   histos[0]->SetMinimum(1e-6);
   histos[0]->SetMarkerStyle(20);
   histos[0]->Draw("lpe");
   for(Int_t i=1; i<tp->NReac(); ++i) {
      histos[i]->SetLineColor(col[(i+1)%14]);
      histos[i]->SetMarkerStyle(20);
      histos[i]->SetMarkerColor(col[(i+1)%14]);
      histos[i]->Draw("samelpe");
   }
   tmg->Draw("C");
   return 0;
}
