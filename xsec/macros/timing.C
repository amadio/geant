#include "TFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"

#include "TPartIndex.h"
#include "TEXsec.h"

#include <string>

using std::string;

void timing(const char *proc="inElastic", const char *part="proton", int elemin=1, int elemax=92)
{
  gSystem->Load("libXsec");

  const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";

  TFile *fx = TFile::Open(fxsec,"r");
  TEXsec *s = (TEXsec *) fx->Get("O");

  int iproc = TPartIndex::I()->ProcIndex(proc);
  int ipart = TPartIndex::I()->PartIndex(part);
  
  float cpu;
  int jele;
  int jpart;
  int jproc;
  float ene;
  int mprod;
  gSystem->Load("libXsec");
  TFile *f = TFile::Open("timing.root");
  f->ls();
  TTree *t = (TTree*) f->Get("G4time");
  t->Print();
  t->SetBranchAddress("cpu",&cpu);
  t->SetBranchAddress("z",&jele);
  t->SetBranchAddress("part",&jpart);
  t->SetBranchAddress("process",&jproc);
  t->SetBranchAddress("energy",&ene);
  t->SetBranchAddress("nprod",&mprod);
  int nentries = t->GetEntries();
  printf("Entries %d\n",nentries);
  TProfile *hh[184];
  for(int ie=0; ie<92; ++ie)
    hh[ie]=hh[92+ie]=0;
  
  for(int ie=0; ie<nentries; ++ie) {
    t->GetEntry(ie);
    if(jele>=elemin && jele<=elemax) {
      if(!hh[jele-1]) {
        string name(TPartIndex::I()->EleSymb(jele));
        string title = string(part)+" "+proc+" on "+name;
        hh[jele-1]=new TProfile((name+"-t").c_str(),title.c_str(),100,
                                log10(TPartIndex::I()->Emin()),
                                log10(TPartIndex::I()->Emax()));
        hh[jele-1]->GetXaxis()->SetTitle("log10(E [GeV])");
        hh[jele-1]->GetYaxis()->SetTitle("t [ms]");
      }
      if(iproc==jproc&&ipart==jpart) hh[jele-1]->Fill(log10(ene),1000*cpu);
    }
  }

  for(int iele=1; iele<=92; ++iele) {
    if(hh[iele-1]) {
      TCanvas *c = new TCanvas();
      hh[iele-1]->Draw();
    }
    if(hh[92+iele-1]) {
      TCanvas *c = new TCanvas();
      hh[92+iele-1]->Draw();
    }
  }
}
