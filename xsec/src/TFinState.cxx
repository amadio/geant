#include "TFinState.h"
#include "TRandom.h"

int TFinState::fVerbose=0;

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
TFinState::TFinState(int nfstates, const int npart[], const float weight[],
                     const float kerma[], const float en[], const char surv[],
                     const int pid[], const float mom[]):
fNFstates(nfstates),
fNsecs(0),
fNMom(0),
fNpart(new int[fNFstates]),
fWeight(new float[fNFstates]),
fKerma(new float[fNFstates]),
fEn(new float[fNFstates]),
fSurv(new char[fNFstates]),
fPID(0),
fMom(0)
{
  memcpy(fNpart,npart,fNFstates*sizeof(int));
  memcpy(fWeight,weight,fNFstates*sizeof(float));
  memcpy(fKerma,kerma,fNFstates*sizeof(float));
  memcpy(fEn,en,fNFstates*sizeof(float));
  memcpy(fSurv,surv,fNFstates*sizeof(char));
  
  fNsecs = 0;
  for(int j=0; j<fNFstates; ++j) {
    fNsecs+=fNpart[j];
    if(j) fWeight[j]+=fWeight[j-1];
  }
  fNMom = 3*fNsecs;

  double wnorm = 1/fWeight[fNFstates-1];
  for(int j=0; j<fNFstates; ++j) fWeight[j]*=wnorm;
  
  fPID = new int[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(int));
  fMom = new float[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(float));
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
    fNpart = new int[fNFstates];
    memcpy(fNpart,right.fNpart,fNFstates*sizeof(int));
    
    delete [] fWeight;
    fWeight = new float[fNFstates];
    memcpy(fWeight,right.fWeight,fNFstates*sizeof(float));
    
    delete [] fKerma;
    fKerma = new float[fNFstates];
    memcpy(fKerma,right.fKerma,fNFstates*sizeof(float));
    
    delete [] fEn;
    fEn = new float[fNFstates];
    memcpy(fEn,right.fEn,fNFstates*sizeof(float));
    
    delete [] fSurv;
    fSurv = new char[fNFstates];
    memcpy(fSurv,right.fSurv,fNFstates*sizeof(char));
    
    delete [] fPID;
    fPID = new int[fNsecs];
    memcpy(fPID,right.fPID,fNsecs*sizeof(int));
    
    delete [] fMom;
    fMom = new float[fNMom];
    memcpy(fMom,right.fMom,fNMom*sizeof(float));
  }
  return *this;
}

//_________________________________________________________________________
bool TFinState::SetFinState(int nfstates, const int npart[], const float weight[],
                              const float kerma[], const float en[], const char surv[],
                              const int pid[], const float mom[])
{
  fNFstates = nfstates;
  
  delete [] fNpart;
  fNpart = new int[fNFstates];
  memcpy(fNpart,npart,fNFstates*sizeof(int));
  
  delete [] fWeight;
  fWeight = new float[fNFstates];
  memcpy(fWeight,weight,fNFstates*sizeof(float));
  
  delete [] fKerma;
  fKerma = new float[fNFstates];
  memcpy(fKerma,kerma,fNFstates*sizeof(float));
  
  delete [] fEn;
  fEn = new float[fNFstates];
  memcpy(fEn,en,fNFstates*sizeof(float));
  
  delete [] fSurv;
  fSurv = new char[fNFstates];
  memcpy(fSurv,surv,fNFstates*sizeof(char));
  
  fNsecs = 0;
  for(int j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];
  fNMom = 3*fNsecs;
  
  delete [] fPID;
  fPID = new int[fNsecs];
  memcpy(fPID,pid,fNsecs*sizeof(int));
  
  delete [] fMom;
  fMom = new float[fNMom];
  memcpy(fMom,mom,fNMom*sizeof(float));

  NormFinSateWeights();
  
  return kTRUE;
}

//_________________________________________________________________________
bool TFinState::SetFinState(const TFinState &right)
{
  fNFstates = right.fNFstates;
  
  delete [] fNpart;
  fNpart = new int[fNFstates];
  memcpy(fNpart,right.fNpart,fNFstates*sizeof(int));
  
  delete [] fWeight;
  fWeight = new float[fNFstates];
  memcpy(fWeight,right.fWeight,fNFstates*sizeof(float));
  
  delete [] fKerma;
  fKerma = new float[fNFstates];
  memcpy(fKerma,right.fKerma,fNFstates*sizeof(float));
  
  delete [] fEn;
  fEn = new float[fNFstates];
  memcpy(fEn,right.fEn,fNFstates*sizeof(float));
  
  delete [] fSurv;
  fSurv = new char[fNFstates];
  memcpy(fSurv,right.fSurv,fNFstates*sizeof(char));
  
  fNsecs = 0;
  for(int j=0; j<fNFstates; ++j) fNsecs+=fNpart[j];
  fNMom = 3*fNsecs;
  
  delete [] fPID;
  fPID = new int[fNsecs];
  memcpy(fPID,right.fPID,fNsecs*sizeof(int));
  
  delete [] fMom;
  fMom = new float[fNMom];
  memcpy(fMom,right.fMom,fNMom*sizeof(float));

  NormFinSateWeights();
  
  return kTRUE;
}


//_________________________________________________________________________
void TFinState::NormFinSateWeights(){
  if(fNsecs)
  {
    for(int j=0; j<fNFstates; ++j) {
      if(j) fWeight[j]+=fWeight[j-1];
    }

    double wnorm = 1/fWeight[fNFstates-1];
    for(int j=0; j<fNFstates; ++j) fWeight[j]*=wnorm; 
  }
} 

//_________________________________________________________________________
bool TFinState::SampleReac(int& npart, float& weight, float& kerma,
                             float &en, const int *&pid, const float *&mom) const
{

  if(!fNFstates) { // ensure that nothing happens 
    npart = 0;
    weight = 0;
    //kerma = 0; //keep current value that is the Ekin of the primary
    //en = kerma;
    pid = 0;
    mom = 0;
    if(kerma<=0.0) { // if it is already stopped 
       en = 0;
       kerma = 0;
       return kFALSE;
    } else { 
      en = -1.; // this case can be checked through if(en<0.)  
      kerma = 0; 
      return kTRUE;
    }
  }

  double eta = gRandom->Rndm();
  int finstat = fNFstates-1;
  for(int i=0; i<fNFstates-1; ++i)
    if(eta<fWeight[i]) {
      finstat = i;
      break;
    }

  int ipoint = 0;
  for(int i=0; i<finstat; ++i) ipoint+=fNpart[i];

  npart = fNpart[finstat];
  weight = fWeight[finstat];
  kerma = fKerma[finstat];
  en = fEn[finstat];
//  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
//  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
  pid = &fPID[ipoint];
  mom = &fMom[3*ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
bool TFinState::SampleReac(int& npart, float& weight, float& kerma,
                    float &en, const int *&pid, const float *&mom, double randn) const
{

  if(!fNFstates) { // ensure that nothing happens 
    npart = 0;
    weight = 0;
    //kerma = 0; //keep current value that is the Ekin of the primary
    //en = kerma;
    pid = 0;
    mom = 0;
    if(kerma<=0.0) { // if it is already stopped 
       en = 0;
       kerma = 0;
       return kFALSE;
    } else { 
      en = -1.; // this case can be checked through if(en<0.)  
      kerma = 0; 
      return kTRUE;
    }
  }

  //double eta = gRandom->Rndm();
  int finstat = fNFstates-1;
  for(int i=0; i<fNFstates-1; ++i)
    if(randn < fWeight[i]) {
      finstat = i;
      break;
    }
  int ipoint = 0;
  for(int i=0; i<finstat; ++i) ipoint+=fNpart[i];

  npart = fNpart[finstat];
  weight = fWeight[finstat];
  kerma = fKerma[finstat];
  en = fEn[finstat];
//  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
//  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
  pid = &fPID[ipoint];
  mom = &fMom[3*ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
bool TFinState::GetReac(int finstat, int& npart, float& weight, float& kerma,
                          float &en, const int *&pid, const float *&mom) const
{
  if(!fNFstates) { // ensure that nothing happens 
    npart = 0;
    weight = 0;
    //kerma = 0; //keep current value that is the Ekin of the primary
    //en = kerma;
    pid = 0;
    mom = 0;
    if(kerma<=0.0) { // if it is already stopped 
       en = 0;
       kerma = 0;
       return kFALSE;
    } else { 
      en = -1.; // this case can be checked through if(en<0.)  
      kerma = 0; 
      return kTRUE;
    }
  } else {
    int ipoint = 0;
    for(int i=0; i<finstat; ++i) ipoint+=fNpart[i];
    npart = fNpart[finstat];
    weight = fWeight[finstat];
    kerma = fKerma[finstat];
    en = fEn[finstat];
    //  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
    //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
    pid = &fPID[ipoint];
    mom = &fMom[3*ipoint];
    return fSurv[finstat];
  }
}


