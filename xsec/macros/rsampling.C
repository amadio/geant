#include <TMath.h>
//#include <TPartIndex.h>
//#include <TEXsec.h>
#include <TH1F.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TText.h>
#include <TMarker.h>

Int_t rsampling(const char *el="O",const char *part="proton",Int_t nrep=100000) 
{
   const EColor col[13] = {kBlack,kGray,kRed,kGreen,kBlue, kMagenta, kCyan,
			   kOrange, kSpring, kTeal, kAzure, kViolet, kPink };

   gSystem->Load("libXsec");
   TFile *f = new TFile("xsec.root");
   f->Get("PartIndex");

   Double_t emin = TPartIndex::I()->Emin();
   Double_t emax = TPartIndex::I()->Emax();

   Int_t nbins = 100;
   Double_t edelta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
   // Sampling of reactions
   Int_t nproc = TPartIndex::I()->NProc();
   Double_t *e = new Double_t[nbins];
   Double_t *x = new Double_t[nbins];
   TH1F **histos = new TH1F*[nproc];

   TMultiGraph *tmg= new TMultiGraph("G5","G5");
   TGraph *gr=0;

   TEXsec * o = (TEXsec*)f->Get(el);
   Double_t en=emin;
   for(Int_t nb=0; nb<nbins; ++nb) {
      e[nb] = TMath::Log10(en);
      en*=edelta;
   }
   char title[200];
   Int_t pindex = TPartIndex::I()->PartIndex(part);
   Double_t ymax = -1;
   Double_t ymin = 1e12;
   for(Int_t i=0; i<nproc-1; ++i) {
      memset(x,0,nbins*sizeof(Double_t));
      if(o->XS(pindex,i,1)>=0) {
	 en=emin;
	 for(Int_t nb=0; nb<nbins; ++nb) {
	    x[nb] = o->XS(pindex, i, en);
	    ymax = x[nb]>ymax?x[nb]:ymax;
	    if(x[nb]) ymin = x[nb]<ymin?x[nb]:ymin;
	    en*=edelta;
	 }
	 gr = new TGraph(nbins,e,x);
	 snprintf(title,199,"Sampling of %s on %s",part,el);
	 gr->SetTitle(title);
	 gr->SetTitle(title);
	 gr->SetLineColor(col[i%13]);
	 tmg->Add(gr);
      }
   }
   ymin = ymin<1e-7?1e-7:ymin;

   char name[10];
   for(Int_t i=0; i<nproc; ++i) {
      snprintf(name,9,"h%.2d",i);
      histos[i]= new TH1F(TPartIndex::I()->ProcName(i),title,nbins,
			  TMath::Log10(emin),TMath::Log10(emax));
   }
   Int_t ibin=nbins;
   en=emin;
   while(ibin--) {
      Double_t hnorm = o->XS(pindex,nproc-1,en)/nrep;
      Double_t len = TMath::Log10(en);
      Int_t irep = nrep;
      while(irep--) {
	 Int_t reac = o->SampleReac(pindex, en);
	 if(reac >= nproc-1) printf("Wrong sampling!! %d\n",reac);
	 histos[reac]->Fill(len,hnorm);
      }
      en*=edelta;
   }
   
   TMarker **m = new TMarker*[nproc];
   TText **t = new TText*[nproc];
   Double_t mx = 0.6;
   Double_t my = 0.88;
   Double_t md = 0.03;
   Bool_t fhist = kFALSE;
   Int_t ndelta = 0;
   for(Int_t i=0; i<nproc; ++i) {
      if(!histos[i]->GetEntries()) continue;
      Int_t colo = col[i%13];
      histos[i]->SetLineColor(colo);
      histos[i]->SetMarkerStyle(20+(ndelta+1)/5);
      histos[i]->SetMarkerColor(colo);
      histos[i]->SetStats(kFALSE);
      if(!fhist) {
	 histos[i]->SetMinimum(ymin);
	 histos[i]->SetMaximum(1.5*ymax);
	 histos[i]->Draw("lpe");
	 fhist = kTRUE;
      }
      else histos[i]->Draw("samelpe");
      m[i] = new TMarker(mx,my-md*ndelta,20+(ndelta+1)/5);
      m[i]->SetNDC();
      m[i]->SetMarkerColor(colo);
      m[i]->Draw();
      t[i] = new TText(mx*1.05,my-md*ndelta,TPartIndex::I()->ProcName(i));
      t[i]->SetTextSize(0.03);
      t[i]->SetNDC();
      t[i]->SetTextAlign(12);
      t[i]->Draw();
      ++ndelta;
   }
   tmg->Draw("C");

   return 0;
}
