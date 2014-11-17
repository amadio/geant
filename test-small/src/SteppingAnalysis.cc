#include <string>

#include "SteppingAnalysis.hh"
//#include "HistogramManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"

#include "G4Material.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eMultipleScattering.hh"

#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4GammaConversion.hh"

#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hIonisation.hh"
#include "G4hPairProduction.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronElasticProcess.hh"

#include "G4MaterialCutsCouple.hh"

#include "G4TrackVector.hh"
#include "TabulatedProcess.hh"

SteppingAnalysis::SteppingAnalysis():
   theHisto(0)

{
  //theHisto->BookHistograms();
}

SteppingAnalysis::~SteppingAnalysis() 
{
//  if(theHisto) delete theHisto;
}

void SteppingAnalysis::DoIt(const G4Step * theStep) 
{
  
  FillCrossSections(theStep);
}

void SteppingAnalysis::FillCrossSections( const G4Step * theStep ) {

  //material check
  //  const G4Material* mat = theStep->GetTrack()->GetMaterial();

#ifdef USE_ROOT

  //fill analysis histograms

  const G4Track* atrack = theStep->GetTrack();
  const G4ParticleDefinition* part = atrack->GetDefinition();

  G4ProcessManager* pm = atrack->GetDefinition()->GetProcessManager();
  G4ProcessVector* pv = pm->GetProcessList();

  int np = pv->length(); 

  G4bool isApp = false;
  G4double meanFreePath = -1.0;
  G4double dedx = -1.0;
  G4double previousStepSize = 0.0;
  G4ForceCondition aCondition;

  const G4MaterialCutsCouple* couple  = atrack->GetMaterialCutsCouple(); 
  G4double  kineticEnergy = atrack->GetKineticEnergy();
  G4float logKineticEnergy = log10(kineticEnergy/GeV);

  //Secondary paticles
  const G4TrackVector* fSecondaryVector = theStep->GetSecondary();
  G4double secondaryKE;

  for(int i = 0; i < np ; ++i) {
    G4VProcess* proc = (*pv)[i];

    G4ProcessType procType = proc->GetProcessType(); 
    const G4String& procName = proc->GetProcessName();
    isApp = proc->IsApplicable(*part);

    if (procType == fElectromagnetic && isApp) {
      //electrons
      if (procName == "eBrem") {
	G4eBremsstrahlung* eproc = (G4eBremsstrahlung*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_eBrem->Fill(logKineticEnergy,1.0/meanFreePath);
	}
	dedx = eproc->GetDEDX(kineticEnergy,couple);
	theHisto->p_dedx_eBrem->Fill(logKineticEnergy,dedx);
	
	// secondaries
	for (unsigned int isec = 0 ; isec < fSecondaryVector->size() ; isec++){
	  G4Track* fSecondaryTrack = (*fSecondaryVector)[isec];
	  secondaryKE = fSecondaryTrack->GetKineticEnergy()/GeV;
	  theHisto->h_secKE_eBrem->Fill(log10(secondaryKE));
	}
      }
      if (procName == "eBremeBrem") {
	TabulatedProcess* eproc = (TabulatedProcess*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_eBrem->Fill(logKineticEnergy,1.0/meanFreePath);
	}
	/*
	dedx = eproc->GetDEDX(kineticEnergy,couple);
	theHisto->p_dedx_eBrem->Fill(logKineticEnergy,dedx);
	
	// secondaries
	for (unsigned int isec = 0 ; isec < fSecondaryVector->size() ; isec++){
	  G4Track* fSecondaryTrack = (*fSecondaryVector)[isec];
	  secondaryKE = fSecondaryTrack->GetKineticEnergy()/GeV;
	  theHisto->h_secKE_eBrem->Fill(log10(secondaryKE));
	}
	*/
      }
      if (procName == "eIoni") {
	G4eIonisation* eproc = (G4eIonisation*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if (meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_eIoni->Fill(logKineticEnergy,1.0/meanFreePath);
	}
	dedx = eproc->GetDEDX(kineticEnergy,couple);
	theHisto->p_dedx_eIoni->Fill(logKineticEnergy,dedx);
      }

      if (procName == "eIonieIoni") {
	TabulatedProcess* eproc = (TabulatedProcess*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if (meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_eIoni->Fill(logKineticEnergy,1.0/meanFreePath);
	}
	//	dedx = eproc->GetDEDX(kineticEnergy,couple);
	//	theHisto->p_dedx_eIoni->Fill(logKineticEnergy,dedx);
      }
      if (procName == "msc") {
	// G4eMultipleScattering* eproc = (G4eMultipleScattering*) proc;
      }

      //photons
      if(procName == "compt") {
	G4ComptonScattering* eproc = (G4ComptonScattering*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_compt->Fill(logKineticEnergy,1.0/meanFreePath);
	}
	
	// secondaries
	for (unsigned int isec = 0 ; isec < fSecondaryVector->size() ; isec++){
	  G4Track* fSecondaryTrack = (*fSecondaryVector)[isec];
	  secondaryKE = fSecondaryTrack->GetKineticEnergy()/GeV;
	  theHisto->h_secKE_compt->Fill(log10(secondaryKE));
	}
      }
      if(procName == "phot") {
	G4PhotoElectricEffect* eproc = (G4PhotoElectricEffect*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_phot->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
      if(procName == "conv") {
	G4GammaConversion* eproc = (G4GammaConversion*) proc;
	meanFreePath = eproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_conv->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
    }
    else if (procType == fHadronic) {
      //hadronic

      if(procName == "hBrems") {
	G4hBremsstrahlung* hproc = (G4hBremsstrahlung*) proc;
	meanFreePath = 
	  hproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_hBrems->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
      if(procName == "hIoni") {
	G4hIonisation* hproc = (G4hIonisation*) proc;
	meanFreePath = 
	  hproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_hIoni->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
      if(procName == "hPairProd") {
	G4hPairProduction* hproc = (G4hPairProduction*) proc;
	meanFreePath = 
	  hproc->MeanFreePath(*atrack);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_hPairProd->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
      if(procName == "hadElastic") {
	G4HadronElasticProcess* hproc = (G4HadronElasticProcess*) proc;
	meanFreePath = 
	  hproc->GetMeanFreePath(*atrack,previousStepSize,&aCondition);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_hadElastic->Fill(logKineticEnergy,1.0/meanFreePath);
	}
      }
      if(procName == "pi-Inelastic") {
	G4PionMinusInelasticProcess* hproc 
	  = (G4PionMinusInelasticProcess*) proc;
	meanFreePath = 
	  hproc->GetMeanFreePath(*atrack,previousStepSize,&aCondition);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_pim_Inelastic->Fill(logKineticEnergy,
					       1./meanFreePath);
	  theHisto->p_nsec_pim_Inelastic->Fill(logKineticEnergy,
					       fSecondaryVector->size());

	  for (unsigned int isec = 0; isec < fSecondaryVector->size();isec++) {
	    G4Track* fSecondaryTrack = (*fSecondaryVector)[isec];
	    secondaryKE = fSecondaryTrack->GetKineticEnergy()/GeV;
	    if(secondaryKE>0.0) {
	      theHisto->p_esec_pim_Inelastic->Fill(logKineticEnergy,
						   secondaryKE);
	      theHisto->h_esec_pim_Inelastic->Fill(log10(secondaryKE/GeV));
	    }
	  }
	}
      }
      if(procName == "pi+Inelastic") {
	G4PionPlusInelasticProcess* hproc = (G4PionPlusInelasticProcess*) proc;
	meanFreePath = 
	  hproc->GetMeanFreePath(*atrack,previousStepSize,&aCondition);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_pip_Inelastic->Fill(logKineticEnergy,
					       1./meanFreePath);
	}
      }
      if(procName == "neutronInelastic") {
	G4NeutronInelasticProcess* hproc = ( G4NeutronInelasticProcess*) proc;
	meanFreePath = 
	  hproc->GetMeanFreePath(*atrack,previousStepSize,&aCondition);
	if(meanFreePath > 0 && meanFreePath < DBL_MAX) {
	  theHisto->p_xsec_N_Inelastic->Fill(logKineticEnergy,1./meanFreePath);
	}
      }
    }
  }
#endif
}

