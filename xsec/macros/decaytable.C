void decaytable(const char *part="proton",Int_t samp=-1)
{
  gSystem->Load("libXsec");
  TFile *fx = new TFile("xsec.root","r");
  TEXsec *s = (TEXsec *) fx->Get("O");
  TFile *ff = new TFile("fstate.root","r");
//  ff->ls();

  Int_t minpart=0;
  Int_t maxpart=TPartIndex::I()->NPart();
  if(TString("*")!=TString(part)) {
     minpart=TPartIndex::I()->PartIndex(part);
     if(minpart<0) {
	printf("Unknown particle\n");
	return;
     }
     maxpart=minpart+1;
  }

  TPDecay *dt = (TPDecay*) ff->Get("DecayTable");
  // Reaction list

  Int_t minsamp = samp;
  Int_t maxsamp = samp;
  if(samp<0 || samp>=dt->NSample()) {
     minsamp = 0;
     maxsamp = dt->NSample();
  }
  
  for(Int_t ipart=minpart; ipart<maxpart; ++ipart) {
     printf("%s\n",TPartIndex::I()->PartName(ipart));
     Int_t npart=0;
     Int_t *pid=0;
     Float_t *mom=0;
     for(Int_t is=minsamp; is<maxsamp; ++is) {
	Bool_t succ = dt->GetDecay(ipart, is, npart, pid, mom);
	if(npart) {
	   Double_t sumpx=0;
	   Double_t sumpy=0;
	   Double_t sumpz=0;
	   Double_t sumen=-TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(ipart))->Mass();
	   Int_t sumch=-TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(ipart))->Charge();
	   for(Int_t ip=0; ip<npart; ++ip) {
	      sumpx+=mom[3*ip  ];
	      sumpy+=mom[3*ip+1];
	      sumpx+=mom[3*ip+2];
	      Double_t dmass=TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(pid[ip]))->Mass();
	      sumen+=TMath::Sqrt(mom[3*ip  ]*mom[3*ip  ]+mom[3*ip+1]*mom[3*ip+1]+
				 mom[3*ip+2]*mom[3*ip+2]+dmass*dmass);
	      sumch+=TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(pid[ip]))->Charge();
	   }
	   if(TMath::Abs(sumpx)+TMath::Abs(sumpy)+TMath::Abs(sumpz)+TMath::Abs(sumen)>1e-5) {
	      printf("#%d np %d: ",is,npart);
	      for(Int_t ip=0; ip<npart; ++ip) printf("%s ",TPartIndex::I()->PartName(pid[ip]));
	      printf("-----------> %g %g %g %g\n",sumpx,sumpy,sumpz,sumen);
	   }
	}
     }
  }
}
