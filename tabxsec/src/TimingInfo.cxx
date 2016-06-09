#include "TFile.h"
#include "TTree.h"
#include "TPartIndex.h"

static TFile *f=0;
static TTree *t=0;

void TimingInfo(float cpu,int z,int pdg,int proc,float en, int np)
{
  static float scpu=0;
  static int sz=0;
  static int spcode=0;
  static int sreac=0;
  static float sen=0;
  static int snp=0;
  static bool first=true;
  
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
    first = false;
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
