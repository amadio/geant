#ifndef GXHistogramManager_H
#define GXHistogramManager_H 1

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TBranch.h"
#endif

class GXHistogramManager 
{
public:

  static GXHistogramManager* Instance();

  GXHistogramManager();
  ~GXHistogramManager();

public:
  void BookHistograms(std::string fileName="test.root");

#ifdef USE_ROOT
  //root histograms
  TFile *theHistFile;

//e bremsstrahlung
  TH1F *h_einc;
  TH1F *h_ntrial;

  TH1F *h_energy[5];
  TH1F *h_angle[5];

  TH1F *g_energy[5];
  TH1F *g_angle[5];

  TH1F *c_energy[5];
  TH1F *c_angle[5];
#endif

private:
  static GXHistogramManager* theInstance;
};

#endif
