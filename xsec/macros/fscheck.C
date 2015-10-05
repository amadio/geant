#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "TPartIndex.h"
#include "TEFstate.h"

using std::string;

enum Error_t {kPart, kProcess};

void usage(Error_t err);

void fscheck(const char *proc="inElastic", const char *part="proton", int elemin=10, int elemax=30)
{
  const char *fxsec = "../../data/xsec_FTFP_BERT_G496p02_1mev.root";
  const char *ffins = "../../data/fstate_FTFP_BERT_G496p02_1mev.root";
  gSystem->Load("libXsec");
  TFile *fx = new TFile(fxsec,"r");
  fx->Get("PartIndex");
  TFile *ff = new TFile(ffins,"r");
//  ff->ls();
  int ireac = TPartIndex::I()->ProcIndex(proc);
  if(ireac<0) {
    printf("Unknown process %s\n",proc);
    usage(kProcess);
    return;
  }
  
  int ipart = TPartIndex::I()->PartIndex(part);
  if(ipart<0) {
    printf("Uknown particle %s\n",part);
    usage(kProcess);
    return;
  }
  
/*
  // Reaction list
  for(int i=0; i<TPartIndex::I()->NProc(); ++i) {
    printf("Proc #%d is %s\n",i,TPartIndex::I()->ProcName(i));
  }
  
  printf("Particles with reactions %d\n",TPartIndex::I()->NPartReac());
  
  for(int i=0; i<TPartIndex::I()->NPartReac();++i) {
    printf("Particle #%d is %s\n",i,TPartIndex::I()->PartName(i));
  }
 */
  
  TEFstate *fs=0;
  for(int i=0; i<92; ++i) {
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(i+1));
//    printf("Element %d is %s\n",i+1,TPartIndex::I()->EleName(i+1));
  }
  
  TProfile *hh[184];
  memset(hh,0,184*sizeof(TProfile*));
  int nh=0;

  char hname_mult[128];
  char hname_pt[128];  
  gStyle->SetOptStat(0);

  for(int iele=elemin-1; iele<elemax; ++iele) {
    hh[iele]=0;
    hh[92+iele]=0;
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(iele+1));
    if(fs) {
      int nsamp = fs->NEFstat();
      string title = string(TPartIndex::I()->PartName(ipart))
      +string(" ")+ string(TPartIndex::I()->ProcName(ireac))+string(" on ")+string(TPartIndex::I()->EleSymb(iele+1));
      sprintf(hname_mult,"%s-mult",TPartIndex::I()->EleSymb(iele+1));
      sprintf(hname_pt,"%s-pt",TPartIndex::I()->EleSymb(iele+1));

      hh[iele] = new TProfile(hname_mult, title.c_str(),
                              100,log10(TPartIndex::I()->Emin()),
                              log10(TPartIndex::I()->Emax()));
      hh[92+iele] = new TProfile(hname_pt, title.c_str(),
                                100,log10(TPartIndex::I()->Emin()),
                                log10(TPartIndex::I()->Emax()));
      hh[iele]->GetXaxis()->SetTitle("Log10(E [GeV])");
      hh[iele]->GetYaxis()->SetTitle("n particles");
      hh[iele]->SetMarkerStyle(20);
      hh[iele]->SetMarkerSize(0.5);
      hh[92+iele]->GetXaxis()->SetTitle("Log10(E [GeV])");
      hh[92+iele]->GetYaxis()->SetTitle("p_{t}");
      hh[92+iele]->SetMarkerStyle(20);
      hh[92+iele]->SetMarkerSize(0.5);
      float kerma=0;
      float weight=0;
      float enr=0;
      int npart=0;
      const int *pid=0;
      const float *mom=0;
      // Test of energy conservation in inelastic
      const double *egrid = TPartIndex::I()->EGrid();
      for(int ien=0; ien<TPartIndex::I()->NEbins(); ++ien) {
        for(int is=0; is<nsamp; ++is) {
          bool isurv = fs->GetReac(ipart, ireac, egrid[ien], is, npart, weight, kerma, enr, pid, mom);
          hh[iele]->Fill(log10(egrid[ien]),npart);
//          printf("Here %d %p\n",npart,mom);
          if(npart)
            hh[92+iele]->Fill(log10(egrid[ien]),sqrt(mom[3*is]*mom[3*is]+mom[3*is+1]*mom[3*is+1]));
	  //          if(is<1000) {
	  //     printf("en %10.3e #%d(%d) ",egrid[ien], npart,isurv);
          //  for(int l=0; l<npart; ++l) if(pid[l]<1000) printf("%s ",TPartIndex::I()->PartName(pid[l]));
	  //    else printf("%d ",pid[l]);
	  //   printf("\n");
	  //  }
        }
      }
    }
  }
  for(int iele=elemin-1; iele<elemax; ++iele) {
    if(hh[iele]) {
      TCanvas *c = new TCanvas();
      hh[iele]->Draw();
    }
    if(hh[92+iele]) {
      TCanvas *c = new TCanvas();
      hh[92+iele]->Draw();
    }
  }
}

void usage(Error_t err)
{
  printf("Usage: finstate.C(const char* process, const char* particle, zmin=1, zmax=92)\n");
  if(err==kProcess) {
    printf("       Available processes: ");
    for(int jp=0; jp<TPartIndex::I()->NProc(); ++jp) {
      printf("%s, ",TPartIndex::I()->ProcName(jp));
      printf("\n");
    }
  } else if(err==kPart) {
    printf("       Particles: ");
    for(int jp=0; jp<TPartIndex::I()->NPartReac(); ++jp) {
      printf("%s, ",TPartIndex::I()->PartName(jp));
      printf("\n");
    }
  }
}
