#include "GUHistogram.h"

namespace vecphys {

GUHistogram::GUHistogram(std::string fileName)
{
  fHistFile = new TFile(fileName.c_str(),"RECREATE");
  BookHistograms();
}

GUHistogram::~GUHistogram()
{
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
}

void GUHistogram::BookHistograms()
{
  fHistFile->mkdir("Benchmark");
  fHistFile->cd("Benchmark");

  ftime   = new TH1F("ftime","ftime",100,0.,0.1);
  fenergy = new TH1F("fenergy","fenergy",100,0.,1000.);
  fangle  = new TH1F("fangle","fangle",100,0.,1.0);
}

} // end namespace vecphys
