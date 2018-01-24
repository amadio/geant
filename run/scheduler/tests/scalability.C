#include <stdio.h>
#include <stdlib.h>
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TROOT.h"

static int irun=0;

void scalability(Int_t option=0, const char *fname="scalability.txt")
{
// Makes scalability plot out of file scalability.txt
// options:
//   0 - scalability
//   1 - run time
//*** run as: root -b -q scalability.C+
   irun++;
   FILE *fp = fopen(fname, "r");
   if (!fp) {
      printf("file scalability.txt not found\n");
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
   while ((read = getline(&line, &len, fp)) != -1) {
      sline = line;
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
   switch (option) {
     case 0:
       hist = new TH1F(Form("scalability%d",irun), "scalability with number of threads", nres, 0.5, nres+0.5);
       hist->GetYaxis()->SetTitle("speed-up");
       break;
     case 1:
       hist = new TH1F(Form("speed%d",irun), "Total run time", nres, 0.5, nres+0.5);
       hist->GetYaxis()->SetTitle("run time [sec]");
       break;
   }    
           
   hist->SetMarkerColor(irun);
   hist->SetMarkerStyle(24+irun);
   hist->SetMarkerSize(1.0);
   hist->SetStats(false);
   hist->GetXaxis()->SetTitle("# threads");
   double error;
   for (int i=0; i<nres; ++i) {
//      printf("%d: %g +/- %g (%g)\n", nth[i],avg[i], sigma[i], avg[0]/avg[i]);
      switch (option) {
        case 0:
          hist->SetBinContent(nth[i],avg[0]/avg[i]);
          error = (avg[0]/avg[i])*sqrt((sigma[i]*sigma[i])/(avg[i]*avg[i])
                                +(sigma[0]*sigma[0])/(avg[0]*avg[0]));
          break;
        case 1:
          hist->SetBinContent(nth[i],avg[i]);
          error = sigma[i];
	  break;
      }  
					
//      printf("error=%g\n", error);				
      hist->SetBinError(nth[i],error);
   }
   TCanvas *c1 = 0;
   Bool_t found = false;
   c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("results");
   if (c1) found = true;
   else c1 = new TCanvas("results", "Results for different number of threads", 700, 800);
   c1->SetGridx();
   c1->SetGridy();
   if (found) hist->Draw("E1SAME");
   else hist->Draw("E1");
   if (option==0 && irun==1) {
     TF1 *lin = new TF1("linear","x",0,nth[nres-1]);
     lin->SetLineColor(kBlue);
     lin->SetLineStyle(kDashed);
     lin->Draw("SAME");
   } 
   switch (option) {
     case 0:
       c1->SaveAs("speedup.gif");
       break;
     case 1:
       c1->SaveAs("time.gif");
   }      
}
