#include "GUHistogram.h"

namespace vecphys {

GUHistogram::GUHistogram(std::string fileName, double maxE)
{
  fHistFile = new TFile(fileName.c_str(),"RECREATE");
  BookHistograms( maxE );
}

GUHistogram::~GUHistogram()
{
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
}

void GUHistogram::BookHistograms( double maxE )
{
  const char* dir= "Benchmark";
  fHistFile->mkdir(dir); // ("Benchmark");
  fHistFile->cd(dir);    // "Benchmark");

  ftime          = new TH1F("ftime",           "Elapsed time",     100, 0.,  0.001);
  fenergyPrimary = new TH1F("Energy_primary",  "Energy Primary",   100, 0.,  1.1 * maxE);  
  fenergyGam     = new TH1F("Energy_gammaOut", "Energy Gamma/Out", 100, 0.,  1.1 * maxE);
  fangleGam      = new TH1F("Angle_gamma",     "Angle  Gamma",     100, -1.0,  1.0);
  fenergyElec    = new TH1F("Energy_elec",     "Energy Electron",  100, 0., maxE);
  fangleElec     = new TH1F("Angle_elec",      "Angle  Electron",  100, -1.0,  1.0);  
}

void GUHistogram::RecordHistos( double EnPrimary,
                                double EnFinalGamma,
                                double angleGamma,    
                                double EnElectron,
                                double angleElectron)
{
   fenergyPrimary->Fill(EnPrimary);
   fenergyGam->Fill(EnFinalGamma);
   fangleGam->Fill(angleGamma);
   fenergyElec->Fill(EnElectron);
   fangleElec->Fill(angleElectron);
}

void GUHistogram::RecordTime( double elapsedTime )
{
    ftime->Fill(elapsedTime);
}
   
} // end namespace vecphys
