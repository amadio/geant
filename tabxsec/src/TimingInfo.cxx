#include "TFile.h"
#include "TTree.h"
#include "TPartIndex.h"

static TFile *f=0;
static TTree *t=0;

void TimingInfo(Float_t cpu,Int_t z,Int_t pdg,Int_t proc,Float_t en, Int_t np)
{
  static Float_t scpu=0;
  static Int_t sz=0;
  static Int_t spcode=0;
  static Int_t sreac=0;
  static Float_t sen=0;
  static Int_t snp=0;
  static Bool_t first=kTRUE;
  
  if(first) {
    f = new TFile("timing.root","recreate");
    f->WriteObject(TPartIndex::I(),"PartIndex");
    t = new TTree("G4time","G4 timings");
    t->Branch("cpu",&scpu,"cpu/F");
    t->Branch("z",&sz,"z/I");
    t->Branch("part",&spcode,"part/I");
    t->Branch("process",&sreac,"process/I");
    t->Branch("energy",&sen,"energy/F");
    t->Branch("nprod",&snp,"nprod/I");
    first = kFALSE;
  }
  //create the file, the Tree and a few branches
  scpu = cpu;
  sz = z;
  spcode = TPartIndex::I()->PartIndex(pdg);
  sreac = TPartIndex::I()->ProcIndex(proc);
  sen = en;
  snp = np;
  t->Fill();
}

void CloseTiming()
{
  if(t)
    t->Write();
  if(f) {
    f->Write();
    f->Close();
  }
}
