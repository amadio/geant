void readHits(bool conc_io){
  R__LOAD_LIBRARY(libGeant_v);
  long nhits = 0;
  TString fname = "hits.root";
  if (conc_io) fname = "hits_output.root";
  TFile *file = TFile::Open(fname);
  if (!file) {printf("file %s not found\n", fname.Data()); return; }
  GeantBlock<MyHit> *fBlock = new GeantBlock<MyHit>();
  fBlock->Initialize(16);
  TTree *fTree = 0;
  if (conc_io) {
    fTree = (TTree*)file->Get("Tree");
    fTree->SetBranchAddress("hitblockoutput", &fBlock);
  } else {
    fTree = (TTree*)file->Get("Hits");
    fTree->SetBranchAddress("hitblocks", &fBlock);
  }    
  const MyHit *hit;
  static int nentries = fTree->GetEntries();
  for (auto entry=0; entry<nentries; ++entry) {
    fTree->GetEntry(entry);
    for (auto ihit=0; ihit<fBlock->Size(); ++ihit) {
      hit = fBlock->At(ihit);
      nhits++;
    }  
  }
  printf("File %s: size = %.2f [MB] compression = %.2f\n", fname.Data(),
        (float)(file->GetSize())/1048576., file->GetCompressionFactor());
  printf("ReadHits read %ld hits\n", nhits);
}
