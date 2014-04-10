#ifndef HistogramManager_H
#define HistogramManager_H 1

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TBranch.h"
#endif

class G4Step;

class HistogramManager 
{
public:

  static HistogramManager* Instance();

  HistogramManager();
  ~HistogramManager();

public:
  void BookHistograms(std::string fileName="TestSmall.root");

#ifdef USE_ROOT
  //root histograms
  TFile *theHistFile;

  TProfile *p_xsec_eBrem;
  TProfile *p_xsec_eIoni;
  TProfile *p_xsec_compt;
  TProfile *p_xsec_phot;
  TProfile *p_xsec_conv;

  TProfile *p_xsec_hBrems;
  TProfile *p_xsec_hIoni;
  TProfile *p_xsec_hadElastic;
  TProfile *p_xsec_hPairProd;
  TProfile *p_xsec_pim_Inelastic;
  TProfile *p_xsec_pip_Inelastic;
  TProfile *p_xsec_N_Inelastic;

  TProfile *p_dedx_eBrem;
  TProfile *p_dedx_eIoni;

  TProfile *p_nsec_pim_Inelastic;
  TProfile *p_esec_pim_Inelastic;
  TH1F *h_esec_pim_Inelastic;

  //Secondaries
  TH1F *h_secKE_eBrem;
  TH1F *h_secKE_compt;

private:
  static HistogramManager* theInstance;

#endif
};

#endif


