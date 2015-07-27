#include "TMath.h"
#include "TPDecay.h"
#include "TFinState.h"
#include "TPartIndex.h"

ClassImp(TPDecay)

//___________________________________________________________________
TPDecay::TPDecay()
: fNSamp(0), fNPart(0), fDecay(0), fCTauPerMass(0)
   //    ,fDecayLambdaTable(0)
{}

//___________________________________________________________________
TPDecay::TPDecay(int nsample, int npart, TFinState *decay)
    : fNSamp(nsample), fNPart(npart), fDecay(decay), fCTauPerMass(0)
//    ,fDecayLambdaTable(0)
{}

//___________________________________________________________________
TPDecay::~TPDecay() {
  //  int npart = TPartIndex::I()->NPart();
  //  for(int i=0; i<npart; ++i)
  //    if(fDecayLambdaTable[i])
  //     delete [] fDecayLambdaTable[i];
  //  delete [] fDecayLambdaTable;
}

//___________________________________________________________________
bool TPDecay::SampleDecay(int pindex, int &npart, const int *&pid, const float *&mom) const {
  float kerma;
  float weight;
  float en;
  return fDecay[pindex].SampleReac(npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
bool TPDecay::GetDecay(int pindex, int ifs, int &npart, const int *&pid, const float *&mom) const {
  float kerma;
  float weight;
  float en;
  return fDecay[pindex].GetReac(ifs, npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
bool TPDecay::HasDecay(int pindex) const {
  if (fDecay[pindex].GetNsecs() == 0)
    return kFALSE;

  return kTRUE;
}

void TPDecay::SetCTauPerMass(double *ctaupermass, int np) {
  if (!fCTauPerMass)
    delete fCTauPerMass;
  fCTauPerMass = new double[np];
  for (int ip = 0; ip < np; ++ip)
    fCTauPerMass[ip] = ctaupermass[ip];
}
