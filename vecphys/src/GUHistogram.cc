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

void GUHistogram::BookHistograms(double maxE)
{
  for(int i = 0 ; i < kNumberPhysicsModel ; ++i) {

    fHistFile->mkdir(GUPhysicsModelName[i]);
    fHistFile->cd(GUPhysicsModelName[i]);

    fTime[i]       = new TH1F("Time",     "Time",     100, 0.,  0.001);
    fEnergyIn[i]   = new TH1F("EnergyIn", "EnergyIn", 100, 0.,  1.1 * maxE);
    fEnergyOut1[i] = new TH1F("EnergyOut1","EnergyOut1", 100, 0.,  1.1 * maxE);
    fEnergyOut2[i] = new TH1F("EnergyOut2","EnergyOut2", 100, 0.,  1.1 * maxE);
    fAngleOut1[i]  = new TH1F("AngleOut1","AngleOut1", 100, -1.,  1.0);
    fAngleOut2[i]  = new TH1F("AngleOut2","AngleOut2", 100, -1.,  1.0);
  }
}

void GUHistogram::RecordHistos(int imodel,
                               double energyIn,
                               double energyOut1,
                               double AngleOut1,
                               double energyOut2,
                               double AngleOut2)
{
  fEnergyIn[imodel]->Fill(energyIn);
  fEnergyOut1[imodel]->Fill(energyOut1);
  fEnergyOut2[imodel]->Fill(energyOut2);
  fAngleOut1[imodel]->Fill(AngleOut1);
  fAngleOut2[imodel]->Fill(AngleOut2);
}

void GUHistogram::RecordTime(int imodel, double elapsedTime)
{
  fTime[imodel]->Fill(elapsedTime);
}

} // end namespace vecphys
