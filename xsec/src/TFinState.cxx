#include "TFinState.h"
#ifdef USE_ROOT
#include "TRandom.h"
#else
#include "base/RNG.h"
using vecgeom::RNG;
#endif

int TFinState::fVerbose = 0;

//_________________________________________________________________________
TFinState::TFinState()
    : fNFstates(0), fNsecs(0), fNMom(0), fWeight(0), fKerma(0), fEn(0), fMom(0), fPID(0), fNpart(0), fSurv(0) {
}

//_________________________________________________________________________
TFinState::TFinState(int nfstates, const int npart[], const float weight[], const float kerma[], const float en[],
                     const char surv[], const int pid[], const float mom[])
    : fNFstates(nfstates), fNsecs(0), fNMom(0), fWeight(new float[fNFstates]), fKerma(new float[fNFstates]),
      fEn(new float[fNFstates]), fMom(0), fPID(0), fNpart(new int[fNFstates]), fSurv(new char[fNFstates]) {
  memcpy(fNpart, npart, fNFstates * sizeof(int));
  memcpy(fWeight, weight, fNFstates * sizeof(float));
  memcpy(fKerma, kerma, fNFstates * sizeof(float));
  memcpy(fEn, en, fNFstates * sizeof(float));
  memcpy(fSurv, surv, fNFstates * sizeof(char));

  fNsecs = 0;
  for (int j = 0; j < fNFstates; ++j) {
    fNsecs += fNpart[j];
    if (j)
      fWeight[j] += fWeight[j - 1];
  }
  fNMom = 3 * fNsecs;

  double wnorm = 1 / fWeight[fNFstates - 1];
  for (int j = 0; j < fNFstates; ++j)
    fWeight[j] *= wnorm;

  fPID = new int[fNsecs];
  memcpy(fPID, pid, fNsecs * sizeof(int));
  fMom = new float[fNMom];
  memcpy(fMom, mom, fNMom * sizeof(float));
}

//_________________________________________________________________________
TFinState::TFinState(const TFinState& other)
    : fNFstates(other.fNFstates), fNsecs(other.fNsecs), fNMom(other.fNMom),
      fWeight(other.fWeight), fKerma(other.fKerma), fEn(other.fEn), fMom(other.fMom),
      fPID(other.fPID), fNpart(other.fNpart), fSurv(other.fSurv) {
}



//_________________________________________________________________________
TFinState::~TFinState() {
  delete[] fNpart;
  delete[] fWeight;
  delete[] fKerma;
  delete[] fEn;
  delete[] fSurv;
  delete[] fPID;
  delete[] fMom;
}

//_________________________________________________________________________
TFinState &TFinState::operator=(const TFinState &right) {
  if (this != &right) {
    fNFstates = right.fNFstates;
    fNsecs = right.fNsecs;
    fNMom = right.fNMom;

    delete[] fNpart;
    fNpart = new int[fNFstates];
    memcpy(fNpart, right.fNpart, fNFstates * sizeof(int));

    delete[] fWeight;
    fWeight = new float[fNFstates];
    memcpy(fWeight, right.fWeight, fNFstates * sizeof(float));

    delete[] fKerma;
    fKerma = new float[fNFstates];
    memcpy(fKerma, right.fKerma, fNFstates * sizeof(float));

    delete[] fEn;
    fEn = new float[fNFstates];
    memcpy(fEn, right.fEn, fNFstates * sizeof(float));

    delete[] fSurv;
    fSurv = new char[fNFstates];
    memcpy(fSurv, right.fSurv, fNFstates * sizeof(char));

    delete[] fPID;
    fPID = new int[fNsecs];
    memcpy(fPID, right.fPID, fNsecs * sizeof(int));

    delete[] fMom;
    fMom = new float[fNMom];
    memcpy(fMom, right.fMom, fNMom * sizeof(float));
  }
  return *this;
}

//_________________________________________________________________________
bool TFinState::SetFinState(int nfstates, const int npart[], const float weight[], const float kerma[],
                            const float en[], const char surv[], const int pid[], const float mom[]) {
  fNFstates = nfstates;

  delete[] fNpart;
  fNpart = new int[fNFstates];
  memcpy(fNpart, npart, fNFstates * sizeof(int));

  delete[] fWeight;
  fWeight = new float[fNFstates];
  memcpy(fWeight, weight, fNFstates * sizeof(float));

  delete[] fKerma;
  fKerma = new float[fNFstates];
  memcpy(fKerma, kerma, fNFstates * sizeof(float));

  delete[] fEn;
  fEn = new float[fNFstates];
  memcpy(fEn, en, fNFstates * sizeof(float));

  delete[] fSurv;
  fSurv = new char[fNFstates];
  memcpy(fSurv, surv, fNFstates * sizeof(char));

  fNsecs = 0;
  for (int j = 0; j < fNFstates; ++j)
    fNsecs += fNpart[j];
  fNMom = 3 * fNsecs;

  delete[] fPID;
  fPID = new int[fNsecs];
  memcpy(fPID, pid, fNsecs * sizeof(int));

  delete[] fMom;
  fMom = new float[fNMom];
  memcpy(fMom, mom, fNMom * sizeof(float));

  NormFinSateWeights();

  return true;
}

//_________________________________________________________________________
bool TFinState::SetFinState(const TFinState &right) {
  fNFstates = right.fNFstates;

  delete[] fNpart;
  fNpart = new int[fNFstates];
  memcpy(fNpart, right.fNpart, fNFstates * sizeof(int));

  delete[] fWeight;
  fWeight = new float[fNFstates];
  memcpy(fWeight, right.fWeight, fNFstates * sizeof(float));

  delete[] fKerma;
  fKerma = new float[fNFstates];
  memcpy(fKerma, right.fKerma, fNFstates * sizeof(float));

  delete[] fEn;
  fEn = new float[fNFstates];
  memcpy(fEn, right.fEn, fNFstates * sizeof(float));

  delete[] fSurv;
  fSurv = new char[fNFstates];
  memcpy(fSurv, right.fSurv, fNFstates * sizeof(char));

  fNsecs = 0;
  for (int j = 0; j < fNFstates; ++j)
    fNsecs += fNpart[j];
  fNMom = 3 * fNsecs;

  delete[] fPID;
  fPID = new int[fNsecs];
  memcpy(fPID, right.fPID, fNsecs * sizeof(int));

  delete[] fMom;
  fMom = new float[fNMom];
  memcpy(fMom, right.fMom, fNMom * sizeof(float));

  NormFinSateWeights();

  return true;
}

//_________________________________________________________________________
void TFinState::NormFinSateWeights() {
  if (fNsecs) {
    for (int j = 0; j < fNFstates; ++j) {
      if (j)
        fWeight[j] += fWeight[j - 1];
    }

    double wnorm = 1 / fWeight[fNFstates - 1];
    for (int j = 0; j < fNFstates; ++j)
      fWeight[j] *= wnorm;
  }
}

//_________________________________________________________________________
bool TFinState::SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid,
                           const float *&mom) const {

  if (!fNFstates) { // ensure that nothing happens
    npart = 0;
    weight = 0;
    // kerma = 0; //keep current value that is the Ekin of the primary
    // en = kerma;
    pid = 0;
    mom = 0;
    if (kerma <= 0.0) { // if it is already stopped
      en = 0;
      kerma = 0;
      return false;
    } else {
      en = -1.; // this case can be checked through if(en<0.)
      kerma = 0;
      return true;
    }
  }

#ifdef USE_ROOT
  double eta = gRandom->Rndm();
#else
  double eta = RNG::Instance().uniform();
#endif
  int finstat = fNFstates - 1;
  for (int i = 0; i < fNFstates - 1; ++i)
    if (eta < fWeight[i]) {
      finstat = i;
      break;
    }

  int ipoint = 0;
  for (int i = 0; i < finstat; ++i)
    ipoint += fNpart[i];

  npart = fNpart[finstat];
  weight = fWeight[finstat];
  kerma = fKerma[finstat];
  en = fEn[finstat];
  //  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
  //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
  pid = &fPID[ipoint];
  mom = &fMom[3 * ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
bool TFinState::SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid, const float *&mom,
                           double randn) const {

  if (!fNFstates) { // ensure that nothing happens
    npart = 0;
    weight = 0;
    // kerma = 0; //keep current value that is the Ekin of the primary
    // en = kerma;
    pid = 0;
    mom = 0;
    if (kerma <= 0.0) { // if it is already stopped
      en = 0;
      kerma = 0;
      return false;
    } else {
      en = -1.; // this case can be checked through if(en<0.)
      kerma = 0;
      return true;
    }
  }

  // double eta = gRandom->Rndm();
  int finstat = fNFstates - 1;
  for (int i = 0; i < fNFstates - 1; ++i)
    if (randn < fWeight[i]) {
      finstat = i;
      break;
    }
  int ipoint = 0;
  for (int i = 0; i < finstat; ++i)
    ipoint += fNpart[i];

  npart = fNpart[finstat];
  weight = fWeight[finstat];
  kerma = fKerma[finstat];
  en = fEn[finstat];
  //  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
  //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
  pid = &fPID[ipoint];
  mom = &fMom[3 * ipoint];
  return fSurv[finstat];
}

//_________________________________________________________________________
bool TFinState::GetReac(int finstat, int &npart, float &weight, float &kerma, float &en, const int *&pid,
                        const float *&mom) const {
  if (!fNFstates) { // ensure that nothing happens
    npart = 0;
    weight = 0;
    // kerma = 0; //keep current value that is the Ekin of the primary
    // en = kerma;
    pid = 0;
    mom = 0;
    if (kerma <= 0.0) { // if it is already stopped
      en = 0;
      kerma = 0;
      return false;
    } else {
      en = -1.; // this case can be checked through if(en<0.)
      kerma = 0;
      return true;
    }
  } else {
    int ipoint = 0;
    for (int i = 0; i < finstat; ++i)
      ipoint += fNpart[i];
    npart = fNpart[finstat];
    weight = fWeight[finstat];
    kerma = fKerma[finstat];
    en = fEn[finstat];
    //  memcpy(pid,&fPID[ipoint],npart*sizeof(int));
    //  memcpy(mom,&fMom[3*ipoint],3*npart*sizeof(float));
    pid = &fPID[ipoint];
    mom = &fMom[3 * ipoint];
    return fSurv[finstat];
  }
}

//_________________________________________________________________________
int TFinState::SizeOf() const {
   size_t size = sizeof(*this);
   size += 3 * fNFstates * sizeof(float);
   size += fNMom * sizeof(float);
   size += fNsecs * sizeof(int);
   size += fNFstates * sizeof(int);
   size += fNFstates * sizeof(char);
   /*
     std::cout << " fNFstates " << fNFstates
	     << " fNsecs " << fNsecs
	     << " fNMom " << fNMom
	     << " size " << (int) size-sizeof(float) << std::endl;
   */
   size -= sizeof(char);  // fStore already takes one char
   size = sizeof(double)*((size-1)/sizeof(double)+1);
   return (int) size;
}

//_________________________________________________________________________
void TFinState::Compact() {
   int size = 0;
   char *start = fStore;
   if(fWeight) {
      size = fNFstates * sizeof(float);
      memcpy(start, fWeight, size);
      fWeight = (float*) start;
      start += size;
   }
   if(fKerma) {
      size = fNFstates * sizeof(float);
      memcpy(start, fKerma, size);
      fKerma = (float*) start;
      start += size;
   }
   if(fEn) {
      size = fNFstates * sizeof(float);
      memcpy(start, fEn, size);
      fEn = (float*) start;
      start += size;
   }
   if(fMom) {
      size = fNMom * sizeof(float);
      memcpy(start, fMom, size);
      fMom = (float*) start;
      start += size;
   }
   if(fPID) {
      size = fNsecs * sizeof(int);
      memcpy(start, fPID, size);
      fPID = (int *) start;
      start += size;
   }
   if(fNpart) {
      size = fNFstates * sizeof(int);
      memcpy(start, fNpart, size);
      fNpart = (int*) start;
      start += size;
   }
   if(fSurv) {
      size = fNFstates * sizeof(char);
      memcpy(start, fSurv, size);
      fSurv = (char*) start;
      start += size;
   }
}

//______________________________________________________________________________
void TFinState::RebuildClass() {
  if(((unsigned long) this) % sizeof(double) != 0) {
      cout << "TFinState::RebuildClass: the class is misaligned" << endl;
      exit(1);
  }
   int size = 0;
   char *start = fStore;
   if(fWeight) {
      //      cout << "Original fWeight " << fWeight << " new pointer " << start << endl;
      fWeight = (float *) start;
      size = fNFstates * sizeof(float);
      start +=size;
   }
   if(fKerma) {
      //      cout << "Original fKerma " << fKerma << " new pointer " << start << endl;
      fKerma = (float *) start;
      size = fNFstates * sizeof(float);
      start +=size;
   }
   if(fEn) {
      //      cout << "Original fEn " << fEn << " new pointer " << start << endl;
      fEn = (float *) start;
      size = fNFstates * sizeof(float);
      start +=size;
   }
   if(fMom) {
      //      cout << "Original fMom " << fMom << " new pointer " << start << endl;
      fMom = (float *) start;
      size = fNMom * sizeof(float);
      start +=size;
   }
   if(fPID) {
      //      cout << "Original fPID " << fPID << " new pointer " << start << endl;
      fPID = (int *) start;
      size = fNsecs * sizeof(int);
      start +=size;
   }
   if(fNpart) {
      //      cout << "Original fNpart " << fNpart << " new pointer " << start << endl;
      fNpart = (int *) start;
      size = fNFstates * sizeof(int);
      start +=size;
   }
   if(fSurv) {
      //      cout << "Original fSurv " << fSurv << " new pointer " << start << endl;
      fSurv = (char *) start;
      size = fNFstates * sizeof(char);
      start +=size;
   }
}
