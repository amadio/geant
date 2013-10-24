#include "TFinState.h"
#include <TRandom.h>

Int_t TFinState::fVerbose=0;

ClassImp(TFinState)

//_________________________________________________________________________
TFinState::TFinState():
fNFstates(0),
fNsecs(0),
fNMom(0),
fPID(0),
fSurv(0),
fNpart(0),
fWeight(0),
fKerma(0),
fMom(0)
{
}

//_________________________________________________________________________
TFinState::TFinState(Int_t nfstates, const Float_t weight[], const Float_t kerma[],
                     const Int_t npart[], const Float_t mom[], const Int_t pid[],
                     const Char_t surv[]):
fNFstates(nfstates),
fNsecs(0),
fNMom(0),
fPID(0),
fSurv(new Char_t[fNFstates]),
fNpart(new Int_t[fNFstates]),
fWeight(new Float_t[fNFstates]),
fKerma(new Float_t[fNFstates]),
fMom(0)
{
  memcpy(fSurv,surv,fNFstates*sizeof(Char_t));
  memcpy(fNpart,npart,fNFstates*sizeof(Int_t));
  memcpy(fWeight,weight,fNFstates*sizeof(Float_t));
  memcpy(fKerma,kerma,fNFstates*sizeof(Float_t));
  
  fNsecs = 0;
  for(Int_t j=0; j<fNFstates; ++j) {
    fNsecs+=fNpart[j];
    if(j) fWeight[j]+=fWeight[j-1];
  }
  fNMom = 3*fNsecs;

  Double_t wnorm = 1/fWeight[fNFstates-1];
  for(Int_t j=0; j<fNFstates; ++j) fWeight[j]*=wnorm;
  
  fMom = new Float_t[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(Float_t));
  fPID = new Int_t[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(Int_t));
}

//_________________________________________________________________________
TFinState::~TFinState()
{
  delete [] fMom;
  delete [] fPID;
  delete [] fSurv;
  delete [] fNpart;
  delete [] fWeight;
  delete [] fKerma;
}

//_________________________________________________________________________
TFinState& TFinState::operator=(const TFinState& right)
{
  if(this != &right) {
    fNFstates = right.fNFstates;
    fNsecs = right.fNsecs;
    fNMom = right.fNMom;
    
    delete [] fPID;
    fPID = new Int_t[fNsecs];
    memcpy(fPID,right.fPID,fNsecs*sizeof(Int_t));
    
    delete [] fSurv;
    fSurv = new Char_t[fNFstates];
    memcpy(fSurv,right.fSurv,fNFstates*sizeof(Char_t));
    
    delete [] fNpart;
    fNpart = new Int_t[fNFstates];
    memcpy(fNpart,right.fNpart,fNFstates*sizeof(Int_t));
    
    delete [] fWeight;
    fWeight = new Float_t[fNFstates];
    memcpy(fWeight,right.fWeight,fNFstates*sizeof(Float_t));
    
    delete [] fKerma;
    fKerma = new Float_t[fNFstates];
    memcpy(fKerma,right.fKerma,fNFstates*sizeof(Float_t));
    
    delete [] fMom;
    fMom = new Float_t[fNMom];
    memcpy(fMom,right.fMom,fNMom*sizeof(Float_t));
  }
  return *this;
}

//_________________________________________________________________________
Bool_t TFinState::SetFinState(Int_t nfstates, const Float_t weight[],
                              const Float_t kerma[], const Int_t npart[],
                              const Float_t mom[], const Int_t pid[], const Char_t surv[])
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
  
  fNsecs = 0;
  for(Int_t j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];
  fNMom = 3*fNsecs;
  
  delete [] fMom;
  fMom = new Float_t[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(Float_t));
  
  delete [] fPID;
  fPID = new Int_t[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(Int_t));
  
  delete [] fSurv;
  fSurv = new Char_t[fNFstates];
  memcpy(fSurv,surv,fNFstates*sizeof(Char_t));
  
  return kTRUE;
}

//_________________________________________________________________________
Bool_t TFinState::SampleReac(Float_t& kerma, Int_t& npart, const Int_t *pid, const Float_t *mom) const
{
  Double_t eta = gRandom->Rndm();
  Int_t finstat = fNFstates-1;
  for(Int_t i=0; i<fNFstates-1; ++i)
    if(eta<fWeight[i]) {
      finstat = i;
      break;
    }
  Int_t ipoint = 0;
  for(Int_t i=0; i<finstat-1; ++i) ipoint+=fNpart[i];
  npart = fNpart[finstat];
  kerma = fKerma[finstat];
//  memcpy(pid,&fPID[ipoint],npart*sizeof(Int_t));
//  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(Float_t));
  pid = &fPID[ipoint];
  mom = &fMom[3*ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
Bool_t TFinState::GetReac(Int_t finstat, Float_t& kerma, Int_t& npart, const Int_t *&pid, const Float_t *mom) const
{
  if(!fNFstates) {
    kerma = 0;
    npart = 0;
    pid = 0;
    mom = 0;
    return kFALSE;
  } else {
    Int_t ipoint = 0;
    for(Int_t i=0; i<finstat-1; ++i) ipoint+=fNpart[i];
    npart = fNpart[finstat];
    kerma = fKerma[finstat];
    //  memcpy(pid,&fPID[ipoint],npart*sizeof(Int_t));
    //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(Float_t));
    pid = &fPID[ipoint];
    mom = &fMom[3*ipoint];
    return fSurv[finstat];
  }
}


