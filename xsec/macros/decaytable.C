#include <string>

#include "TFile.h"
#include "TSystem.h"

#include "TEXsec.h"
#include "TPdecay.h"

using std::string;

void decaytable(const char *part="pi-",int samp=-1)
{
  gSystem->Load("libXsec");

  const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
  const char *ffins = "../../data/fstate_FTFP_BERT_G496p02_1mev.root";

  TFile *fx = new TFile(fxsec,"r");
  TEXsec *s = (TEXsec *) fx->Get("O");
  TFile *ff = new TFile(ffins,"r");
//  ff->ls();

  int minpart=0;
  int maxpart=TPartIndex::I()->NPart();
  if(string("*")!=string(part)) {
     minpart=TPartIndex::I()->PartIndex(part);
     if(minpart<0) {
	printf("Unknown particle\n");
	return;
     }
     maxpart=minpart+1;
  }

  TPDecay *dt = (TPDecay*) ff->Get("DecayTable");
  // Reaction list

  int minsamp = samp;
  int maxsamp = samp;
  if(samp<0 || samp>=dt->NSample()) {
     minsamp = 0;
     maxsamp = dt->NSample();
  }
  
  for(int ipart=minpart; ipart<maxpart; ++ipart) {
     printf("%s\n",TPartIndex::I()->PartName(ipart));
     int npart=0;
     const int *pid=0;
     const float *mom=0;
     for(int is=minsamp; is<maxsamp; ++is) {
	bool succ = dt->GetDecay(ipart, is, npart, pid, mom);
	if(npart) {
	   double sumpx=0;
	   double sumpy=0;
	   double sumpz=0;
	   double sumen=-TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(ipart))->Mass();
	   int sumch=-TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(ipart))->Charge();
	   for(int ip=0; ip<npart; ++ip) {
	      sumpx+=mom[3*ip  ];
	      sumpy+=mom[3*ip+1];
	      sumpx+=mom[3*ip+2];
	      double dmass=TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(pid[ip]))->Mass();
	      sumen+=sqrt(mom[3*ip  ]*mom[3*ip  ]+mom[3*ip+1]*mom[3*ip+1]+
				 mom[3*ip+2]*mom[3*ip+2]+dmass*dmass);
	      sumch+=TDatabasePDG::Instance()->GetParticle(TPartIndex::I()->PDG(pid[ip]))->Charge();
	   }
	   if(fabs(sumpx)+fabs(sumpy)+fabs(sumpz)+fabs(sumen)>1e-5) {
	      printf("#%d np %d: ",is,npart);
	      for(int ip=0; ip<npart; ++ip) printf("%s ",TPartIndex::I()->PartName(pid[ip]));
	      printf("-----------> %g %g %g %g\n",sumpx,sumpy,sumpz,sumen);
	   }
	}
     }
  }
}
