void finstate(const char *proc="inElastic", const char *part="proton", Int_t elemin=1, Int_t elemax=92)
{
  gSystem->Load("libXsec");
  TFile *fx = new TFile("xsec.root","r");
  TEXsec *s = (TEXsec *) fx->Get("O");
  TFile *ff = new TFile("fstate.root","r");
//  ff->ls();
  iproc = TPartIndex::I()->ProcIndex(proc);
  ipart = TPartIndex::I()->PartIndex(part);

  
  // Reaction list
  for(Int_t i=0; i<TPartIndex::I()->NProc(); ++i) {
    printf("Proc #%d is %s\n",i,TPartIndex::I()->ProcName(i));
  }
  
  printf("Particles with reactions %d\n",TPartIndex::I()->NPartReac());
  
  for(Int_t i=0; i<TPartIndex::I()->NPartReac();++i) {
    printf("Particle #%d is %s\n",i,TPartIndex::I()->PartName(i));
  }
  
  TEFstate *fs=0;
  for(Int_t i=0; i<92; ++i) {
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(i+1));
//    printf("Element %d is %s\n",i+1,TPartIndex::I()->EleName(i+1));
  }
  
  TProfile *hh[184];
  Int_t nh=0;
  
  for(Int_t iele=elemin; iele<=elemax; ++iele) {
    hh[iele-1]=0;
    hh[92+iele-1]=0;
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(iele));
    if(fs) {
      Int_t nsamp = fs->NEFstat();
      TString title = TString(TPartIndex::I()->PartName(ipart))
      +TString(" ")+ TString(TPartIndex::I()->ProcName(iproc))+TString(" on ")+TString(TPartIndex::I()->EleSymb(iele));
      TString name = TString(TPartIndex::I()->EleSymb(iele));
      hh[iele-1] = new TProfile(name+"-mult", title,
                              100,TMath::Log10(TPartIndex::I()->Emin()),
                              TMath::Log10(TPartIndex::I()->Emax()));
      hh[92+iele-1] = new TProfile(name+"-pt", title,
                                100,TMath::Log10(TPartIndex::I()->Emin()),
                                TMath::Log10(TPartIndex::I()->Emax()));
      hh[iele-1]->GetXaxis()->SetTitle("Log10(GeV)");
      hh[iele-1]->GetYaxis()->SetTitle("n particles");
      hh[92+iele-1]->GetXaxis()->SetTitle("Log10(GeV)");
      hh[92+iele-1]->GetYaxis()->SetTitle("pt");
      Float_t kerma=0;
      Int_t npart=0;
      Int_t *pid=0;
      Float_t *mom=0;
      // Test of energy conservation in inelastic
      Double_t *egrid = TPartIndex::I()->EGrid();
      for(Int_t ien=0; ien<TPartIndex::I()->NEbins(); ++ien) {
        for(Int_t is=0; is<nsamp; ++is) {
          Int_t isurv = fs->GetReac(ipart, egrid[ien], iproc, is, kerma, npart, pid, mom);
          hh[iele-1]->Fill(TMath::Log10(egrid[ien]),npart);
//          printf("Here %d %p\n",npart,mom);
          if(npart)
            hh[92+iele-1]->Fill(TMath::Log10(egrid[ien]),TMath::Sqrt(mom[3*is]*mom[3*is]+mom[3*is+1]*mom[3*is+1]));
          if(is<1000) {
//            printf("en %10.3e #%d(%d) ",egrid[ien], npart,isurv);
//            for(Int_t l=0; l<npart; ++l) if(pid[l]<1000) printf("%s ",TPartIndex::I()->PartName(pid[l]));
//            else printf("%d ",pid[l]);
//            printf("\n");
          }
        }
      }
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
