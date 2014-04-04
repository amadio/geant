#include "TFinState.h"
#include "TRandom.h"

Int_t TFinState::fVerbose=0;

ClassImp(TFinState)

//_________________________________________________________________________
TFinState::TFinState():
fNFstates(0),
fNsecs(0),
fNMom(0),
fNpart(0),
fWeight(0),
fKerma(0),
fEn(0),
fSurv(0),
fPID(0),
fMom(0)
{
}

//_________________________________________________________________________
TFinState::TFinState(Int_t nfstates, const Int_t npart[], const Float_t weight[],
                     const Float_t kerma[], const Float_t en[], const Char_t surv[],
                     const Int_t pid[], const Float_t mom[]):
fNFstates(nfstates),
fNsecs(0),
fNMom(0),
fNpart(new Int_t[fNFstates]),
fWeight(new Float_t[fNFstates]),
fKerma(new Float_t[fNFstates]),
fEn(new Float_t[fNFstates]),
fSurv(new Char_t[fNFstates]),
fPID(0),
fMom(0)
{
  memcpy(fNpart,npart,fNFstates*sizeof(Int_t));
  memcpy(fWeight,weight,fNFstates*sizeof(Float_t));
  memcpy(fKerma,kerma,fNFstates*sizeof(Float_t));
  memcpy(fEn,en,fNFstates*sizeof(Float_t));
  memcpy(fSurv,surv,fNFstates*sizeof(Char_t));
  
  fNsecs = 0;
  for(Int_t j=0; j<fNFstates; ++j) {
    fNsecs+=fNpart[j];
    if(j) fWeight[j]+=fWeight[j-1];
  }
  fNMom = 3*fNsecs;

  Double_t wnorm = 1/fWeight[fNFstates-1];
  for(Int_t j=0; j<fNFstates; ++j) fWeight[j]*=wnorm;
  
  fPID = new Int_t[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(Int_t));
  fMom = new Float_t[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(Float_t));
}

//_________________________________________________________________________
TFinState::~TFinState()
{
  delete [] fNpart;
  delete [] fWeight;
  delete [] fKerma;
  delete [] fEn;
  delete [] fSurv;
  delete [] fPID;
  delete [] fMom;
}

//_________________________________________________________________________
TFinState& TFinState::operator=(const TFinState& right)
{
  if(this != &right) {
    fNFstates = right.fNFstates;
    fNsecs = right.fNsecs;
    fNMom = right.fNMom;
    
    delete [] fNpart;
    fNpart = new Int_t[fNFstates];
    memcpy(fNpart,right.fNpart,fNFstates*sizeof(Int_t));
    
    delete [] fWeight;
    fWeight = new Float_t[fNFstates];
    memcpy(fWeight,right.fWeight,fNFstates*sizeof(Float_t));
    
    delete [] fKerma;
    fKerma = new Float_t[fNFstates];
    memcpy(fKerma,right.fKerma,fNFstates*sizeof(Float_t));
    
    delete [] fEn;
    fEn = new Float_t[fNFstates];
    memcpy(fEn,right.fEn,fNFstates*sizeof(Float_t));
    
    delete [] fSurv;
    fSurv = new Char_t[fNFstates];
    memcpy(fSurv,right.fSurv,fNFstates*sizeof(Char_t));
    
    delete [] fPID;
    fPID = new Int_t[fNsecs];
    memcpy(fPID,right.fPID,fNsecs*sizeof(Int_t));
    
    delete [] fMom;
    fMom = new Float_t[fNMom];
    memcpy(fMom,right.fMom,fNMom*sizeof(Float_t));
  }
  return *this;
}

//_________________________________________________________________________
Bool_t TFinState::SetFinState(Int_t nfstates, const Int_t npart[], const Float_t weight[],
                              const Float_t kerma[], const Float_t en[], const Char_t surv[],
                              const Int_t pid[], const Float_t mom[])
{
  fNFstates = nfstates;
  
  delete [] fNpart;
  fNpart = new Int_t[fNFstates];
  memcpy(fNpart,npart,fNFstates*sizeof(Int_t));
  
  delete [] fWeight;
  fWeight = new Float_t[fNFstates];
  memcpy(fWeight,weight,fNFstates*sizeof(Float_t));
  
  delete [] fKerma;
  fKerma = new Float_t[fNFstates];
  memcpy(fKerma,kerma,fNFstates*sizeof(Float_t));
  
  delete [] fEn;
  fEn = new Float_t[fNFstates];
  memcpy(fEn,en,fNFstates*sizeof(Float_t));
  
  delete [] fSurv;
  fSurv = new Char_t[fNFstates];
  memcpy(fSurv,surv,fNFstates*sizeof(Char_t));
  
  fNsecs = 0;
  for(Int_t j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];
  fNMom = 3*fNsecs;
  
  delete [] fPID;
  fPID = new Int_t[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(Int_t));
  
  delete [] fMom;
  fMom = new Float_t[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(Float_t));

  NormFinSateWeights();
  
  return kTRUE;
}

//_________________________________________________________________________
void TFinState::NormFinSateWeights(){
  for(Int_t j=0; j<fNFstates; ++j) {
    if(j) fWeight[j]+=fWeight[j-1];
  }

  Double_t wnorm = 1/fWeight[fNFstates-1];
  for(Int_t j=0; j<fNFstates; ++j) fWeight[j]*=wnorm; 
} 

//_________________________________________________________________________
Bool_t TFinState::SampleReac(Int_t& npart, Float_t& weight, Float_t& kerma,
                             Float_t &en, const Int_t *&pid, const Float_t *&mom) const
{
  Double_t eta = gRandom->Rndm();
  Int_t finstat = fNFstates-1;
  for(Int_t i=0; i<fNFstates-1; ++i)
    if(eta<fWeight[i]) {
      finstat = i;
      break;
    }
  Int_t ipoint = 0;
  for(Int_t i=0; i<finstat; ++i) ipoint+=fNpart[i];

  npart = fNpart[finstat];
  weight = fWeight[finstat];
  kerma = fKerma[finstat];
  en = fEn[finstat];
//  memcpy(pid,&fPID[ipoint],npart*sizeof(Int_t));
//  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(Float_t));
  pid = &fPID[ipoint];
  mom = &fMom[3*ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
Bool_t TFinState::GetReac(Int_t finstat, Int_t& npart, Float_t& weight, Float_t& kerma,
                          Float_t &en, const Int_t *&pid, const Float_t *&mom) const
{
  if(!fNFstates) {
    npart = 0;
    weight = 0;
    kerma = 0;
    en = 0;
    pid = 0;
    mom = 0;
    return kFALSE;
  } else {
    Int_t ipoint = 0;
    for(Int_t i=0; i<finstat; ++i) ipoint+=fNpart[i];
    npart = fNpart[finstat];
    weight = fWeight[finstat];
    kerma = fKerma[finstat];
    en = fEn[finstat];
    //  memcpy(pid,&fPID[ipoint],npart*sizeof(Int_t));
    //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(Float_t));
    pid = &fPID[ipoint];
    mom = &fMom[3*ipoint];
    return fSurv[finstat];
  }
}


