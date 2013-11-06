void timing(const char *proc="inElastic", const char *part="proton", Int_t elemin=1, Int_t elemax=92)
{
  gSystem->Load("libXsec");
  TFile *fx = new TFile("xsec.root","r");
  TEXsec *s = (TEXsec *) fx->Get("O");

  Int_t iproc = TPartIndex::I()->ProcIndex(proc);
  Int_t ipart = TPartIndex::I()->PartIndex(part);
  
  Float_t cpu;
  Int_t jele;
  Int_t jpart;
  Int_t jproc;
  Float_t ene;
  Int_t mprod;
  gSystem->Load("libXsec");
  TFile *f = new TFile("timing.root");
  f->ls();
  TTree *t = (TTree*) f->Get("G4time");
  t->Print();
  t->SetBranchAddress("cpu",&cpu);
  t->SetBranchAddress("z",&jele);
  t->SetBranchAddress("part",&jpart);
  t->SetBranchAddress("process",&jproc);
  t->SetBranchAddress("energy",&ene);
  t->SetBranchAddress("nprod",&mprod);
  Int_t nentries = t->GetEntries();
  printf("Entries %d\n",nentries);
  TProfile *hh[184];
  for(Int_t ie=0; ie<92; ++ie)
    hh[ie]=hh[92+ie]=0;
  
  for(Int_t ie=0; ie<nentries; ++ie) {
    t->GetEntry(ie);
    if(jele>=elemin && jele<=elemax) {
      if(!hh[jele-1]) {
        TString name(TPartIndex::I()->EleSymb(jele));
        TString title = TString(part)+" "+proc+" on "+name;
        hh[jele-1]=new TProfile(name+"-t",title,100,
                                TMath::Log10(TPartIndex::I()->Emin()),
                                TMath::Log10(TPartIndex::I()->Emax()));
        hh[jele-1]->GetXaxis()->SetTitle("Log10(GeV)");
        hh[jele-1]->GetYaxis()->SetTitle("t(ms)");
      }
      if(iproc==jproc&&ipart==jpart) hh[jele-1]->Fill(TMath::Log10(ene),1000*cpu);
    }
  }

  for(Int_t iele=1; iele<=92; ++iele) {
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
