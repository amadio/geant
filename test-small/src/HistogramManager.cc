//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "HistogramManager.hh"
#include "TabulatedPhysicsNameSpace.hh"

HistogramManager* HistogramManager::theInstance = 0;

HistogramManager* HistogramManager::Instance()
{
  if (theInstance == 0) theInstance = new HistogramManager();
  return theInstance;
}

HistogramManager::HistogramManager() 
{
}

HistogramManager::~HistogramManager() 
{
#ifdef USE_ROOT
  theHistFile->cd();
  theHistFile->Write();
  theHistFile->Close();
#endif
}

void HistogramManager::BookHistograms(std::string fileName) 
{
#ifdef USE_ROOT

  const int nbin = TP::Nbin;
  const int emin = TP::Emin;
  const int emax = TP::Emax;

  theHistFile = new TFile(fileName.c_str(),"RECREATE");

  theHistFile->mkdir("MeanFreePath");
  theHistFile->cd("MeanFreePath");

  p_xsec_eBrem = new TProfile("p_xsec_eBrem","p_xsec_eBrem",nbin,emin,emax);
  p_xsec_eIoni = new TProfile("p_xsec_eIoni","p_xsec_eIoni",nbin,emin,emax);
  p_xsec_compt = new TProfile("p_xsec_compt","p_xsec_compt",nbin,emin,emax);
  p_xsec_phot  = new TProfile("p_xsec_phot" ,"p_xsec_phot" ,nbin,emin,emax);
  p_xsec_conv  = new TProfile("p_xsec_conv" ,"p_xsec_conv" ,nbin,emin,emax);

  p_xsec_hBrems = new TProfile("p_xsec_hBrems","p_xsec_hBrems",nbin,emin,emax);
  p_xsec_hIoni = new TProfile("p_xsec_hIoni","p_xsec_hIoni",nbin,emin,emax);
  p_xsec_hadElastic = new TProfile("p_xsec_hadElastic",
				   "p_xsec_hadElastic",nbin,emin,emax);
  p_xsec_hPairProd = new TProfile("p_xsec_hPairProd",
				  "p_xsec_hPairProd",nbin,emin,emax);
  p_xsec_pim_Inelastic = new TProfile("p_xsec_pim_Inelastic",
				      "p_xsec_pim_Inelastic" ,nbin,emin,emax);

  p_nsec_pim_Inelastic = new TProfile("p_nsec_pim_Inelastic",
				      "p_nsec_pim_Inelastic" ,nbin,emin,emax);
  p_esec_pim_Inelastic = new TProfile("p_esec_pim_Inelastic",
				      "p_esec_pim_Inelastic" ,nbin,emin,emax);
  h_esec_pim_Inelastic = new TH1F("h_esec_pim_Inelastic",
				  "h_esec_pim_Inelastic" ,nbin,-10,emax);

  p_xsec_pip_Inelastic = new TProfile("p_xsec_pip_Inelastic",
				      "p_xsec_pip_Inelastic" ,nbin,emin,emax);
  p_xsec_N_Inelastic   = new TProfile("p_xsec_N_Inelastic",
				      "p_xsec_N_Inelastic" ,nbin,emin,emax);

  p_dedx_eBrem = new TProfile("p_dedx_eBrem","p_dedx_eBrem",nbin,emin,emax);
  p_dedx_eIoni = new TProfile("p_dedx_eIoni","p_dedx_eIoni",nbin,emin,emax);

  theHistFile->mkdir("Secondary");
  theHistFile->cd("Secondary");

  h_secKE_eBrem = new TH1F("h_secKE_eBrem","h_secKE_eBrem" ,140,-10,4.0);
  h_secKE_compt = new TH1F("h_secKE_compt","h_secKE_compt" ,140,-10,4.0);

#endif
}
