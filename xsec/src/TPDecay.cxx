#include "TMath.h"
#include "TPDecay.h"
#include "TFinState.h"
#include "TPartIndex.h"
//#include <iostream>

ClassImp(TPDecay)

//___________________________________________________________________
TPDecay::TPDecay():
    fNSamp(0),
    fNPart(0),
    fDecay(0),
    fCTauPerMass(0)
//    ,fDecayLambdaTable(0)
{
}

//___________________________________________________________________
TPDecay::TPDecay(Int_t nsample, Int_t npart, TFinState *decay):
    fNSamp(nsample),
    fNPart(npart),
    fDecay(decay),
    fCTauPerMass(0)
//    ,fDecayLambdaTable(0)
{
}

//___________________________________________________________________
TPDecay::~TPDecay() {
//  Int_t npart = TPartIndex::I()->NPart();
//  for(Int_t i=0; i<npart; ++i) 
//    if(fDecayLambdaTable[i]) 
//     delete [] fDecayLambdaTable[i];
//  delete [] fDecayLambdaTable;
}

//___________________________________________________________________
Bool_t TPDecay::SampleDecay(Int_t pindex, Int_t &npart,
                            const Int_t *&pid, const Float_t *&mom) const
{
  Float_t kerma;
  Float_t weight;
  Float_t en;
  return fDecay[pindex].SampleReac(npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
Bool_t TPDecay::GetDecay(Int_t pindex, Int_t ifs, Int_t &npart,
                         const Int_t *&pid, const Float_t *&mom) const
{
  Float_t kerma;
  Float_t weight;
  Float_t en;
  return fDecay[pindex].GetReac(ifs,npart, weight, kerma, en, pid, mom);
}

//___________________________________________________________________
Bool_t TPDecay::HasDecay(Int_t pindex) const {
  if(fDecay[pindex].GetNsecs()==0)
    return kFALSE;

  return kTRUE;
}

void  TPDecay::SetCTauPerMass(Double_t *ctaupermass, Int_t np){
  if(!fCTauPerMass)
    delete fCTauPerMass;
  fCTauPerMass = new Double_t[np];
  for(Int_t ip=0; ip<np; ++ip)
     fCTauPerMass[ip] = ctaupermass[ip];
}

