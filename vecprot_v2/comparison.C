#include <stdio.h>
#include <stdlib.h>
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"

static int irun=0;

void comparison(const char *fname="file.txt")
{
// Makes scalability plot out of file scalability.txt
// options:
//   0 - scalability
//   1 - run time
//*** run as: root -b -q scalability.C+
   irun++;
   FILE *fp = fopen(fname, "r");
   if (!fp) {
      printf("file %s not found\n", fname);
      return;
   }
   
   char * line = NULL;
   int nres = -1;
   double time[100];
   int nth[100];
   double avg[100];
   double sigma[100];
   TString sline;
   size_t len = 0;
   ssize_t read;
   int nthreads = 0;
   int ntot = 0;
   double sum = 0;
   double sumsq = 0;
   double mean = 0;
   TString hname;
   while ((read = getline(&line, &len, fp)) != -1) {
      sline = line;
      if (sline.Contains("#")) {
         hname = sline;
         hname.ReplaceAll("#","");
         continue;
      }   
      if (sline.Contains("NTHREADS")) {
         sline.ReplaceAll("NTHREADS ","");
         if (ntot>0) {
            nth[nres] = nthreads;
            avg[nres] = sum/ntot;   
	         for (int i=0; i<ntot; ++i) {
	            sumsq += (time[i]-avg[nres])*(time[i]-avg[nres]);
            }
            sigma[nres] = sqrt(sumsq/ntot);	    
         }
         nthreads = sline.Atoi();
         nres++;
         ntot = 0;
         sum = 0;
         sumsq = 0;
         continue;
      }
      time[ntot] = sline.Atof();
      sum += time[ntot];
      ntot++;
   }
   nth[nres] = nthreads;
   avg[nres] = sum/ntot;
   for (int i=0; i<ntot; ++i) {
      sumsq += (time[i]-avg[nres])*(time[i]-avg[nres]);
   }
   sigma[nres] = sqrt(sumsq/ntot);
   nres++;
   delete line;

   TH1F *hist = 0;
   hist = new TH1F(Form("plot%d",irun), hname, nres, 0.5, nth[nres-1]+0.5);
   hist->GetYaxis()->SetTitle(hname);
   hist->GetYaxis()->SetRangeUser(4.E5, 1.E8);       
   hist->SetMarkerColor(irun);
   hist->SetFillColor(irun);
//   hist->SetFillStyle(3003+irun);
   hist->SetMarkerStyle(24+irun);
   hist->SetMarkerSize(1.0);
   hist->SetBarWidth(0.2);
   hist->SetBarOffset(0.2*irun);
   hist->SetStats(kFALSE);
   hist->GetXaxis()->SetTitle("# threads");
   double error;
   for (int i=0; i<nres; ++i) {
      hist->SetBinContent(nth[i],avg[i]);
      hist->SetBinError(nth[i],sigma[i]);
   }
   TCanvas *c1 = 0;
   Bool_t found = false;
   c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("results");
   if (c1) found = true;
   else c1 = new TCanvas("results", "Results for different number of threads", 700, 800);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetLogy();
   if (found) {
     hist->Draw("E1BAR SAME");
   }  
   else hist->Draw("E1BAR");
   c1->SaveAs("comparison.gif");
}
