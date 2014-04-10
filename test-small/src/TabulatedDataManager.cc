//
// S.Y. Jun & J. Apostolakis, April 2014
//
#include "TabulatedDataManager.hh"

#include "TSystem.h"
#include "TFile.h"
#include <TPartIndex.h>
#include <TEXsec.h>
#include <TEFstate.h>

TabulatedDataManager::TabulatedDataManager() :
  mxsec(0), 
  fstrack_h(0)
{
  //tabulated physics data
  char* fxsecFileName = getenv("VP_DATA_XSEC");
  char* ffstaFileName = getenv("VP_DATA_FSTA");

  if (fxsecFileName) fxsec = new TFile(fxsecFileName,"r");
  else               fxsec = new TFile("xsec_FTFP.root","r");
  if (ffstaFileName) ffsta = new TFile(ffstaFileName,"r");
  else               ffsta = new TFile("fstate_FTFP.root","r");
}

TabulatedDataManager::~TabulatedDataManager() {
  free(mxsec);
  free(fstrack_h);
}

void TabulatedDataManager::PrepareTable(const char* pnam, const char* reac) 
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
  Int_t nsamp = 10;
  
  mxsec = (double*) malloc((max_element+1)*nbins*sizeof(double));
  fstrack_h = (GXTrack *) malloc (nbins*nsamp*max_element*sizeof(GXTrack));
  
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
      mxsec[ib] += xs*weight[ie];
      mxsec[ib+(ie+1)*nbins] = xs*weight[ie];
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
        if(mom) {
          fstrack_h[dindex].px = mom[3*is];
          fstrack_h[dindex].py = mom[3*is+1];
          fstrack_h[dindex].pz = mom[3*is+2];
          fstrack_h[dindex].E  = enr;
        }
      }
    }
  }
}
