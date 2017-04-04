#include "GUHistogram.h"
#include "base/PhysicalConstants.h"

namespace vecphys {

GUHistogram::GUHistogram(std::string fileName, double maxE)
{
  fHistFile = new TFile(fileName.c_str(), "RECREATE");
  BookHistograms(maxE);
}

GUHistogram::~GUHistogram()
{
  fHistFile->cd();
  fHistFile->Write();
  fHistFile->Close();
}
    
void GUHistogram::BookHistograms(double minE, double maxE)
{
    for (int i = 0; i < kNumberPhysicsModel; ++i)
    {
        fHistFile->mkdir(GUPhysicsModelName[i]);
        fHistFile->cd(GUPhysicsModelName[i]);
        
        fTime[i] = new TH1F("Time", "Time", 100, 0., 0.001);
        fEnergyIn[i] = new TH1F("EnergyIn", "EnergyIn", 100, 0.9*minE, 1.1*maxE);
        fEnergyOut1[i] = new TH1F("EnergyOut1", "EnergyOut1", 100, 0.9*minE, 1.1*maxE);
        fEnergyOut2[i] = new TH1F("EnergyOut2", "EnergyOut2", 100, 0.9*minE, 1.1*maxE);
        fAngleOut1[i] = new TH1F("AngleOut1", "AngleOut1", 100, -1., 1.0);
        fAngleOut2[i] = new TH1F("AngleOut2", "AngleOut2", 100, -1., 1.0);
    }
        
    for (int i = 0; i < kNumberPhysicsProcess; ++i) {
    
        fHistFile->mkdir(GUPhysicsProcessName[i]);
        fHistFile->cd(GUPhysicsProcessName[i]);
        
        fProcEnergy[i] = new TH1F("ProcEnergy", "ProcEnergy", 100, 0.9*minE, 1.1*maxE);
        fNint[i] = new TH1F("Nint", "Nint", 100, 0., 10.);
        if (i == kPhotonProcess) {
            fStep[i] = new TH1F("Step", "Step", 100, 0., 1.);
            fLambda[i] = new TH1F("Lambda", "Lambda", 100, 0., 1.);
        } else {
            fStep[i] = new TH1F("Step", "Step", 100, 0., 100.);
            fLambda[i] = new TH1F("Lambda", "Lambda", 100, 0., 100.);
        }
    }
}

void GUHistogram::BookHistograms(double maxE)
{
  for (int i = 0; i < kNumberPhysicsModel; ++i) {

    fHistFile->mkdir(GUPhysicsModelName[i]);
    fHistFile->cd(GUPhysicsModelName[i]);

    fTime[i] = new TH1F("Time", "Time", 100, 0., 0.001);
    if(i<2){
        double minE= maxE / (1 + 2.0 * maxE * inv_electron_mass_c2);
        fEnergyIn[i] = new TH1F("EnergyIn", "EnergyIn", 100, 0.9*minE, 1.1*maxE);
        fEnergyOut1[i] = new TH1F("EnergyOut1", "EnergyOut1", 100, 0.9*minE, 1.1*maxE);
        fEnergyOut2[i] = new TH1F("EnergyOut2", "EnergyOut2", 100, 0.9*minE, 1.1*maxE);
        
    }
    else {
        fEnergyIn[i] = new TH1F("EnergyIn", "EnergyIn", 100, 0., 1.1 * maxE);
        fEnergyOut1[i] = new TH1F("EnergyOut1", "EnergyOut1", 100, 0., 1.1 * maxE);
        fEnergyOut2[i] = new TH1F("EnergyOut2", "EnergyOut2", 100, 0., 1.1 * maxE);
    }
    fAngleOut1[i] = new TH1F("AngleOut1", "AngleOut1", 100, -1., 1.0);
    fAngleOut2[i] = new TH1F("AngleOut2", "AngleOut2", 100, -1., 1.0);
  }

  for (int i = 0; i < kNumberPhysicsProcess; ++i) {

    fHistFile->mkdir(GUPhysicsProcessName[i]);
    fHistFile->cd(GUPhysicsProcessName[i]);

    fProcEnergy[i] = new TH1F("ProcEnergy", "ProcEnergy", 100, 0., 1.1 * maxE);
    fNint[i] = new TH1F("Nint", "Nint", 100, 0., 10.);
    if (i == kPhotonProcess) {
      fStep[i] = new TH1F("Step", "Step", 100, 0., 1.);
      fLambda[i] = new TH1F("Lambda", "Lambda", 100, 0., 1.);
    } else {
      fStep[i] = new TH1F("Step", "Step", 100, 0., 100.);
      fLambda[i] = new TH1F("Lambda", "Lambda", 100, 0., 100.);
    }
  }
}

void GUHistogram::RecordHistos(int imodel, double energyIn, double energyOut1, double AngleOut1, double energyOut2,
                               double AngleOut2)
{
  fEnergyIn[imodel]->Fill(energyIn);
  fEnergyOut1[imodel]->Fill(energyOut1);
  fEnergyOut2[imodel]->Fill(energyOut2);
  fAngleOut1[imodel]->Fill(AngleOut1);
  fAngleOut2[imodel]->Fill(AngleOut2);
}

void GUHistogram::RecordHistosProc(int iprocess, double energy, double nint, double step, double lambda)
{
  fProcEnergy[iprocess]->Fill(energy);
  fNint[iprocess]->Fill(nint);
  fStep[iprocess]->Fill(step);
  fLambda[iprocess]->Fill(lambda);
}

void GUHistogram::RecordTime(int imodel, double elapsedTime) { fTime[imodel]->Fill(elapsedTime); }

} // end namespace vecphys
