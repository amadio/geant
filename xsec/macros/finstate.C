#include "TSystem.h"
#include "TFile.h"
#include "TProfile.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TEXsec.h"
#include "TEFstate.h"

void finstate(const char *proc="inElastic", const char *part="proton", Int_t elemin=1, Int_t elemax=92)
{
  gSystem->Load("libXsec");
  TFile *fx = new TFile("xsec.root","r");
  fx->Get("PartIndex");
  TFile *ff = new TFile("fstate.root","r");
//  ff->ls();
  Int_t iproc = TPartIndex::I()->ProcIndex(proc);
  Int_t ipart = TPartIndex::I()->PartIndex(part);

  
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
  memset(hh,0,184*sizeof(TProfile*));
  Int_t nh=0;
  
  for(Int_t iele=elemin-1; iele<elemax; ++iele) {
    hh[iele]=0;
    hh[92+iele]=0;
    fs = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(iele+1));
    if(fs) {
      Int_t nsamp = fs->NEFstat();
      TString title = TString(TPartIndex::I()->PartName(ipart))
      +TString(" ")+ TString(TPartIndex::I()->ProcName(iproc))+TString(" on ")+TString(TPartIndex::I()->EleSymb(iele+1));
      TString name = TString(TPartIndex::I()->EleSymb(iele+1));
      hh[iele] = new TProfile(name+"-mult", title,
                              100,TMath::Log10(TPartIndex::I()->Emin()),
                              TMath::Log10(TPartIndex::I()->Emax()));
      hh[92+iele] = new TProfile(name+"-pt", title,
                                100,TMath::Log10(TPartIndex::I()->Emin()),
                                TMath::Log10(TPartIndex::I()->Emax()));
      hh[iele]->GetXaxis()->SetTitle("Log10(GeV)");
      hh[iele]->GetYaxis()->SetTitle("n particles");
      hh[92+iele]->GetXaxis()->SetTitle("Log10(GeV)");
      hh[92+iele]->GetYaxis()->SetTitle("pt");
      Float_t kerma=0;
      Int_t npart=0;
      const Int_t *pid=0;
      const Float_t *mom=0;
      // Test of energy conservation in inelastic
      const Double_t *egrid = TPartIndex::I()->EGrid();
      for(Int_t ien=0; ien<TPartIndex::I()->NEbins(); ++ien) {
        for(Int_t is=0; is<nsamp; ++is) {
          Int_t isurv = fs->GetReac(ipart, egrid[ien], iproc, is, kerma, npart, pid, mom);
          hh[iele]->Fill(TMath::Log10(egrid[ien]),npart);
//          printf("Here %d %p\n",npart,mom);
          if(npart)
            hh[92+iele]->Fill(TMath::Log10(egrid[ien]),TMath::Sqrt(mom[3*is]*mom[3*is]+mom[3*is+1]*mom[3*is+1]));
          if(is<1000) {
            printf("en %10.3e #%d(%d) ",egrid[ien], npart,isurv);
            for(Int_t l=0; l<npart; ++l) if(pid[l]<1000) printf("%s ",TPartIndex::I()->PartName(pid[l]));
            else printf("%d ",pid[l]);
            printf("\n");
          }
        }
      }
    }
  }
  for(Int_t iele=elemin-1; iele<elemax; ++iele) {
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
