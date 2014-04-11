//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "TabulatedProcess.hh"
#include "TabulatedDataManager.hh"
#include "TabulatedPhysicsNameSpace.hh"

#include "G4Track.hh"
#include "G4VParticleChange.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessType.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4GPILSelection.hh"
#include "G4SystemOfUnits.hh"

#include "TSystem.h"
#include <TPartIndex.h>
#include <TEXsec.h>
#include <TEFstate.h>

TFile* TabulatedProcess::fxsec = 0;
TFile* TabulatedProcess::ffsta = 0;

TabulatedProcess::TabulatedProcess(G4String processName, 
				   G4ProcessType processType) :
  G4WrapperProcess(processName,processType),
  particleChange(0),
  fXsecTable(0),
  fFinalState(0)
{

  //tabulated physics data
  char* fxsecFileName = getenv("VP_DATA_XSEC");
  char* ffstaFileName = getenv("VP_DATA_FSTA");

  if (fxsecFileName) {
    fxsec = new TFile(fxsecFileName,"r");
    std::cout << "Using " << fxsecFileName << std::endl;
  }
  else {
    fxsec = new TFile("xsec_FTFP.root","r");
  }

  if (ffstaFileName) {
    ffsta = new TFile(ffstaFileName,"r");
    std::cout << "Using " << ffstaFileName << std::endl;
  }
  else {
    ffsta = new TFile("fstate_FTFP.root","r");
  }

  //  dataManager = new TabulatedDataManager();

  std::cout  << "processName " << processName << std::endl;
  if(processName == "eIoni") {
    PrepareTable("e-","Ionisation"); //process name from VP
  }
  else if (processName == "eBrem") {
    PrepareTable("e-","Brehms"); //process name from VP
  }
}

TabulatedProcess::~TabulatedProcess() 
{
  fxsec->Close();
  ffsta->Close();
  free(fXsecTable);
  free(fFinalState);
}

G4double
TabulatedProcess::PostStepGetPhysicalInteractionLength(const G4Track& track,
						   G4double previousStepSize,
						   G4ForceCondition* condition)
{
  // Call to ensure that the process is correctly initialised - temporary
  pRegProcess->PostStepGetPhysicalInteractionLength(track,
						    previousStepSize,
						    condition);

  if ( (previousStepSize <=0.0) || (theNumberOfInteractionLengthLeft<=0.0)) {
    // beggining of tracking (or just after DoIt of this process)
    ResetNumberOfInteractionLengthLeft();
  } else if ( previousStepSize > 0.0) {
    // subtract NumberOfInteractionLengthLeft 
    SubtractNumberOfInteractionLengthLeft(previousStepSize);
  } else {
    // zero step
    //  DO NOTHING
  }

  // condition is set to "Not Forced"
  *condition = NotForced;

  // get mean free path
  currentInteractionLength = MeanFreePath(track);

  G4double value;
  if (currentInteractionLength <DBL_MAX) {
    value = theNumberOfInteractionLengthLeft * currentInteractionLength;
  } else {
    value = DBL_MAX;
  }
  return value;
}

G4double TabulatedProcess::MeanFreePath(const G4Track& track)
{
  // DefineMaterial(track.GetMaterialCutsCouple());

  G4double logKE = log10(track.GetKineticEnergy()/GeV);
  G4double deltaE = 1.0*(TP::Emax-TP::Emin)/TP::Nbin;

  G4int ibin = int((logKE-TP::Emin)/deltaE);
  if (ibin < 0) ibin = 0;
  if (ibin > TP::Nbin-1) ibin = TP::Nbin -2;

  //put an interpolation within a bin here vl+(vh-vl)*(x-xl)/(xh-xl)
  G4double yL = fXsecTable[ibin];  
  G4double yH = fXsecTable[ibin+1];  

  G4double dx = (logKE-TP::Emin-ibin)/deltaE;  
  G4double preStepLambda = yL + (yH-yL)*dx;
  
  G4double x = DBL_MAX;
  if(0.0 < preStepLambda) { x = 1.0/preStepLambda; }
  return x;
}

G4VParticleChange* 
TabulatedProcess::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  // process PostStepDoIt for the original process
  particleChange = pRegProcess->PostStepDoIt(track, step);

  return particleChange;
}

void TabulatedProcess::Print(const G4Step& step) 
{
  G4StepPoint* postStepPoint = step.GetPostStepPoint();
  std::cout << "TabulatedProcess ProcName, E, dE, x, StepLength Nsec " 
	    << postStepPoint->GetProcessDefinedStep()->GetProcessName() << " " 
	    << postStepPoint->GetKineticEnergy()/MeV  << " " 
	    << step.GetTotalEnergyDeposit()/MeV  << " " 
	    << postStepPoint->GetPosition() << " "  
	    << step.GetStepLength() << " " 
	    << particleChange->GetNumberOfSecondaries() << std::endl;
}

void TabulatedProcess::PrepareTable(const char* pnam, const char* reac) 
{
  //no checks for particle and reaction for now

  TPartIndex* tp = (TPartIndex*)fxsec->Get("PartIndex");
  const Int_t nbins = TPartIndex::I()->NEbins();
  Int_t ipart = TPartIndex::I()->PartIndex(pnam);
  
  //setup for material (PbWO4) - the part should be generalized with G4Material
  const Int_t max_element = 3;
  const Int_t element[max_element] = {82,74,8};
  const Int_t natom[max_element]   = {1,1,4};
  const Float_t watom[max_element] = {207.2,183.84,16};
  
  Float_t tweight = 0.0;
  Float_t weight[max_element]; 
  for(Int_t ie = 0; ie < max_element ; ++ie) tweight += watom[ie]*natom[ie];

  //number of final states sampled by VP
  Int_t nsamp = 5;
  
  fXsecTable = (double*) malloc((max_element+1)*nbins*sizeof(double));
  fFinalState = (GXTrack *) malloc (nbins*nsamp*max_element*sizeof(GXTrack));
  
  for(Int_t ie = 0; ie < max_element ; ++ie) {
    
    const char* mat = TPartIndex::I()->EleSymb(element[ie]);
    
    //cross section per element with the weight in [mm^-1]
    weight[ie] = 6.02486*8.28*natom[ie]/tweight; 
    TEXsec *mate = (TEXsec*)fxsec->Get(mat);
    
    Double_t emin = mate->Emin();
    Double_t emax = mate->Emax();
    Double_t delta = exp(log(emax/emin)/(nbins-1));
    //xsec per element and the composite material
    Int_t rcode = tp->ProcIndex(reac);
    Double_t en=emin;
    for(Int_t ib=0; ib < nbins; ++ib) {
      Float_t xs = mate->XS(ipart,rcode,en);
      fXsecTable[ib] += xs*weight[ie];
      fXsecTable[ib+(ie+1)*nbins] = xs*weight[ie];
      en*=delta;
    }
    TEFstate *fss = (TEFstate *) 
      ffsta->Get(TPartIndex::I()->EleSymb(element[ie]));

    Int_t ireac = TPartIndex::I()->ProcIndex(reac);
    size_t dindex = 0;

    for(Int_t ien=0; ien < TPartIndex::I()->NEbins(); ++ien) {
      for(Int_t is=0; is<nsamp; ++is) {
        dindex = is + nsamp*(ien + nbins+nbins*ie);       
        
        Float_t kerma=0;
        Float_t w=0;
        Float_t enr=0;
        Int_t npart=0;
        const Int_t *pid=0;
        const Float_t *mom=0;
        
        // Test of energy conservation in inelastic
        const Double_t *egrid = TPartIndex::I()->EGrid();
        
        Bool_t isurv = fss->GetReac(ipart, ireac, egrid[ien], is, npart, 
                                    w, kerma, enr, pid, mom);
        if(npart>0) {
          fFinalState[dindex].px = mom[3*is];
          fFinalState[dindex].py = mom[3*is+1];
          fFinalState[dindex].pz = mom[3*is+2];
          fFinalState[dindex].E  = enr;
        }
      }
    }
  }
}
