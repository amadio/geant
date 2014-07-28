#include <stdio.h>
#include <stdlib.h>
#include "TString.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"

void scalability()
{
// Makes scalability plot out of file scalability.txt
//*** run as: root -b -q scalability.C+
   FILE *fp = fopen("scalability.txt", "r");
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

   TH1F *hist = new TH1F("scalability", "scalability with number of threads",
   nres, 0.5, nres+0.5);
   hist->SetMarkerColor(kRed);
   hist->SetMarkerStyle(25);
   hist->SetMarkerSize(0.9);
   hist->SetStats(kFALSE);
   hist->GetXaxis()->SetTitle("# threads");
   hist->GetYaxis()->SetTitle("speed-up");
   double error;
   for (int i=0; i<nres; ++i) {
//      printf("%d: %g +/- %g (%g)\n", nth[i],avg[i], sigma[i], avg[0]/avg[i]);
      hist->SetBinContent(nth[i],avg[0]/avg[i]);
      error = (avg[0]/avg[i])*sqrt((sigma[i]*sigma[i])/(avg[i]*avg[i])
                                +(sigma[0]*sigma[0])/(avg[0]*avg[0]));
//      printf("error=%g\n", error);				
      hist->SetBinError(nth[i],error);
   }   
   TCanvas *c1 = new TCanvas("Speed-up", "Speed-up against nthreads", 700, 800);
   c1->SetGridx();
   c1->SetGridy();
   hist->Draw();
   TF1 *lin = new TF1("linear","x",0,nth[nres-1]);
   lin->SetLineColor(kBlue);
   lin->SetLineStyle(kDashed);
   lin->Draw("SAME");
   c1->SaveAs("speedup.gif");
}
