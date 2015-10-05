#include <TPartIndex.h>
#include <TEXsec.h>
#include <TH1F.h>
#include <TFile.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TText.h>
#include <TMarker.h>

int rsampling(const char *el="O",const char *part="proton",int nrep=100000) 
{
   const EColor col[13] = {kBlack,kGray,kRed,kGreen,kBlue, kMagenta, kCyan,
			   kOrange, kSpring, kTeal, kAzure, kViolet, kPink };

   gSystem->Load("libXsec");

   const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";

   TFile *f = TFile::Open(fxsec,"r");
   f->Get("PartIndex");

   double emin = TPartIndex::I()->Emin();
   double emax = TPartIndex::I()->Emax();

   int nbins = 100;
   double edelta = exp(log(emax/emin)/(nbins-1));
   // Sampling of reactions
   int nproc = TPartIndex::I()->NProc();
   double *e = new double[nbins];
   double *x = new double[nbins];
   TH1F **histos = new TH1F*[nproc];

   TMultiGraph *tmg= new TMultiGraph("GV","GV");
   TGraph *gr=0;

   TEXsec * o = (TEXsec*)f->Get(el);
   double en=emin;
   for(int nb=0; nb<nbins; ++nb) {
      e[nb] = log10(en);
      en*=edelta;
   }
   char title[200];
   int pindex = TPartIndex::I()->PartIndex(part);
   double ymax = -1;
   double ymin = 1e12;
   for(int i=0; i<nproc-1; ++i) {
      memset(x,0,nbins*sizeof(double));
      if(o->XS(pindex,i,1)>=0) {
	 en=emin;
	 for(int nb=0; nb<nbins; ++nb) {
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
   for(int i=0; i<nproc; ++i) {
      snprintf(name,9,"h%.2d",i);
      histos[i]= new TH1F(TPartIndex::I()->ProcName(i),title,nbins,
			  log10(emin),log10(emax));
   }
   int ibin=nbins;
   en=emin;
   while(ibin--) {
      double hnorm = o->XS(pindex,nproc-1,en)/nrep;
      double len = log10(en);
      int irep = nrep;
      while(irep--) {
	 int reac = o->SampleReac(pindex, en);
	 if(reac >= nproc-1) printf("Wrong sampling!! %d\n",reac);
	 histos[reac]->Fill(len,hnorm);
      }
      en*=edelta;
   }
   
   TCanvas *c1 = new TCanvas();
   c1->SetLogy();

   TMarker **m = new TMarker*[nproc];
   TText **t = new TText*[nproc];
   double mx = 0.6;
   double my = 0.88;
   double md = 0.03;
   bool fhist = kFALSE;
   int ndelta = 0;
   for(int i=0; i<nproc; ++i) {
      if(!histos[i]->GetEntries()) continue;
      int colo = col[i%13];
      histos[i]->SetLineColor(colo);
      histos[i]->SetMarkerStyle(20+(ndelta+1)/5);
      histos[i]->SetMarkerColor(colo);
      histos[i]->SetStats(kFALSE);
      if(!fhist) {
	 histos[i]->SetMinimum(ymin);
	 histos[i]->SetMaximum(1.5*ymax);
	 histos[i]->Draw("lpe");
	 histos[i]->GetXaxis()->SetTitle("Log10(E [GeV])");
	 histos[i]->GetYaxis()->SetTitle("#sigma [barn]");
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
