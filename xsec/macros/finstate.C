void finstate() 
{
   gSystem->Load("libXsec");
   TFile *fx = new TFile("xsec.root","r");
   TEXsec *s = (TEXsec *) fx->Get("O");
   TFile *ff = new TFile("fstate.root","r");
   ff->ls();
  
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
     printf("Element %d is %s\n",i+1,TPartIndex::I()->EleName(i+1));
  }
  
  for(Int_t iele=0; iele<92; ++iele) {
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(iele));
    if(fs) {
      Int_t iproc = 4; // inelastic
      Int_t ipart = 32; // proton
      Float_t kerma=0;
      Int_t npart=0;
      Int_t *pid=0;
      Float_t *mom=0;
      // Test of energy conservation in inelastic
      Double_t *egrid = TPartIndex::I()->EGrid();
      
      for(Int_t ien=0; ien<TPartIndex::I()->NEbins(); ++ien) {
        Int_t isurv = fs->GetReac(ipart, egrid[ien], iproc, 0, kerma, npart, pid, mom);
        printf("en %f #%d(%d) ",egrid[ien], npart,isurv);
        for(Int_t l=0; l<npart; ++l) if(pid[l]<1000) printf("%s ",TPartIndex::I()->PartName(pid[l]));
        else printf("%d ",pid[l]);
        printf("\n");
      }
    }
  }
}
