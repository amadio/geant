#include <string>
#include "GXHistogramManager.h"
#include "stdio.h"

GXHistogramManager* GXHistogramManager::theInstance = 0;

GXHistogramManager* GXHistogramManager::Instance()
{
  if (theInstance == 0) theInstance = new GXHistogramManager();
  return theInstance;
}

GXHistogramManager::GXHistogramManager() 
{
  ;
}

GXHistogramManager::~GXHistogramManager() 
{
#ifdef USE_ROOT
  theHistFile->cd();
  theHistFile->Write();
  theHistFile->Close();
#endif
}

void GXHistogramManager::BookHistograms(std::string fileName) 
{
#ifdef USE_ROOT

  printf("Booking Histograms\n");
  const int nbin = 100;
  const int emin = 1.0;
  const int emax = 1001.;

  theHistFile = new TFile(fileName.c_str(),"RECREATE");

  //Bremsstrahlung

  theHistFile->mkdir("Validation");
  theHistFile->cd("Validation");

  h_einc = new TH1F("h_einc","h_einc" ,nbin,emin,emax);
  h_ntrial = new TH1F("h_ntrial","h_trial" ,20,0,20.);

  char hname[128];

  const G4int nModels = 5;
  char* model[nModels] = {"pdf","pdf2","alias","alias2","g4"};

  for(int i = 0 ; i < nModels ; ++i) {
    sprintf(hname,"h_energy_%s",model[i]);
    h_energy[i] = new TH1F(hname,hname,nbin,emin,emax);
    sprintf(hname,"h_angle_%s",model[i]);
    h_angle[i] = new TH1F(hname,hname,nbin,0.0,1.0);
  }

  theHistFile->mkdir("GPU");
  theHistFile->cd("GPU");

  for(int i = 0 ; i < nModels ; ++i) {
    sprintf(hname,"h_energy_%s",model[i]);
    g_energy[i] = new TH1F(hname,hname,nbin,emin,emax);
    sprintf(hname,"h_angle_%s",model[i]);
    g_angle[i] = new TH1F(hname,hname,nbin,0.0,1.0);
  }

  theHistFile->mkdir("CPU");
  theHistFile->cd("CPU");

  for(int i = 0 ; i < nModels ; ++i) {
    sprintf(hname,"c_energy_%s",model[i]);
    c_energy[i] = new TH1F(hname,hname,nbin,emin,emax);
    sprintf(hname,"c_angle_%s",model[i]);
    c_angle[i] = new TH1F(hname,hname,nbin,0.0,1.0);
  }
#endif
}
