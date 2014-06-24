#include "TAxis.h"
#include "TCanvas.h"
#include "TEFstate.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLine.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include <TObjArray.h>
#include <TObjString.h>
#include <TPFstate.h>
#include <TPartIndex.h>
#include <TROOT.h>
#include <TString.h>
#include <TText.h>

ClassImp(TEFstate)


TEFstate* TEFstate::fElements[NELEM]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t TEFstate::fNLdElems=0; 


//___________________________________________________________________
TEFstate::TEFstate():
   fEle(0),
   fDens(0),
   fAtcm3(0),
   fEmin(0),
   fEmax(0),
   fNEbins(0),
   fEilDelta(0),
   fEGrid(TPartIndex::I()->EGrid()),
   fNEFstat(0),
   fNRpart(0),
   fPFstate(0)
{
}

//___________________________________________________________________
TEFstate::TEFstate(Int_t z, Int_t a, Float_t dens):
   fEle(z*10000+a*10),
   fDens(dens),
   fAtcm3(fDens*TMath::Na()*1e-24/TPartIndex::I()->WEle(z)),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fNEbins(TPartIndex::I()->NEbins()),
   fEilDelta(TPartIndex::I()->EilDelta()),
   fEGrid(TPartIndex::I()->EGrid()),
   fNEFstat(0),
   fNRpart(TPartIndex::I()->NPartReac()),
   fPFstate(new TPFstate[fNRpart])
{
}

//___________________________________________________________________
TEFstate::~TEFstate() {
   delete [] fPFstate;
}

//___________________________________________________________________
Bool_t TEFstate::AddPart(Int_t kpart, Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[])
{
   return fPFstate[kpart].SetPart(pdg,nfstat,nreac,dict);
}

//___________________________________________________________________
Bool_t TEFstate::AddPart(Int_t kpart, Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[], TFinState vecfs[])
{
  if(!fNEFstat) fNEFstat = nfstat;
  else if(fNEFstat != nfstat) {
    ::Fatal("AddPart","Number of final sample states changed during run from %d to %d",fNEFstat,nfstat);
  }
  return fPFstate[kpart].SetPart(pdg,nfstat,nreac,dict,vecfs);
}

//___________________________________________________________________
Bool_t TEFstate::AddPartFS(Int_t kpart, Int_t ibin, Int_t reac, const Int_t npart[], const Float_t weight[],
                           const Float_t kerma[], const Float_t en[], const Char_t surv[], const Int_t pid[],
                           const Float_t mom[])
{
  return fPFstate[kpart].SetFinState(ibin, reac, npart, weight, kerma, en, surv, pid, mom);
}

//_____________________________________________________________________________
void TEFstate::SetRestCaptFstate(Int_t kpart, const TFinState &fstate){
	fPFstate[kpart].SetRestCaptFstate(fstate);
}

//______________________________________________________________________________
Bool_t TEFstate::HasRestCapture(Int_t partindex)
{
  if( partindex<TPartIndex::I()->NPartReac() ) 
    return fPFstate[partindex].HasRestCaptFstat();
  return kFALSE; 
}

//______________________________________________________________________________
Bool_t TEFstate::SampleRestCaptFstate(Int_t kpart,Int_t& npart, Float_t& weight,
                            Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const
{
   if(kpart<TPartIndex::I()->NPartReac()){	
   	return fPFstate[kpart].SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom);
   } else {
    kerma=0;
    npart=0;
    pid=0;
    mom=0;
    return kFALSE;
   }	
}

//______________________________________________________________________________
Bool_t TEFstate::SampleRestCaptFstate(Int_t kpart,Int_t& npart, Float_t& weight,
                            Float_t& kerma, Float_t &enr, const Int_t *&pid, 
                            const Float_t *&mom, Double_t randn) const
{
   if(kpart<TPartIndex::I()->NPartReac()){	
   	return fPFstate[kpart].SampleRestCaptFstate(npart, weight, kerma, enr, pid, mom, randn);
   } else {
    kerma=0;
    npart=0;
    pid=0;
    mom=0;
    return kFALSE;
   }	
}

//___________________________________________________________________
Bool_t TEFstate::SampleReac(Int_t pindex, Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                            Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const
{
  return fPFstate[pindex].SampleReac(preac, en, npart, weight, kerma, enr, pid, mom);
}

//___________________________________________________________________
Bool_t TEFstate::SampleReac(Int_t pindex, Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                            Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom,
                            Double_t randn1, Double_t randn2) const
{
  return fPFstate[pindex].SampleReac(preac, en, npart, weight, kerma, enr, pid, mom, randn1, randn2);
}

//___________________________________________________________________
Bool_t TEFstate::GetReac(Int_t pindex, Int_t preac, Float_t en, Int_t ifs, Int_t& npart, Float_t& weight,
                         Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const
{
  return fPFstate[pindex].GetReac(preac, en, ifs, npart, weight, kerma, enr, pid, mom);
}


//___________________________________________________________________
TEFstate *TEFstate::GetElement(Int_t z, Int_t a, TFile* f) {
   //   printf("Getting Element %d %d %d\n",z,a,fNLdElems);
   Int_t ecode = z*10000+a*10;
   for(Int_t el=0; el<fNLdElems; ++el) 
      if(ecode == fElements[el]->Ele()) 
	 return fElements[el];

   // Element not found in memory, getting it from file
   TFile *ff=gFile;
   if(f) ff=f;
   if(!ff) ::Fatal("TEFstate::GetElement","No file open!");
   fElements[fNLdElems] = (TEFstate *) ff->Get(TPartIndex::I()->EleSymb(z));
   if(!fElements[fNLdElems]) {
      ::Fatal("GetElement","Element z %d a %d not found",z,a);
      return 0; // just to make the compiler happy
   } else {
      // We loaded the element, however we have to see whether
      // the energy grid is the right one
      // NO, don't need to check. It will be loaded from xsec.root
//      if(FloatDiff(TPartIndex::I()->Emin(),fElements[fNLdElems]->Emin(),1e-7) ||
//	 FloatDiff(TPartIndex::I()->Emax(),fElements[fNLdElems]->Emax(),1e-7) ||
//	 TPartIndex::I()->NEbins() != fElements[fNLdElems]->NEbins())
	 // we have to resize the energy grid of the element
//	 fElements[fNLdElems]->Resample();
      return fElements[fNLdElems++];
   }
}

//___________________________________________________________________
Bool_t TEFstate::Prune()
{
   for(Int_t ip=0; ip<fNRpart; ++ip)
      fPFstate[ip].Prune();
   return kTRUE;
}

//___________________________________________________________________
Bool_t TEFstate::Resample()
{
   for(Int_t ip=0; ip<fNRpart; ++ip)
      fPFstate[ip].Resample();
   fEmin = TPartIndex::I()->Emin();
   fEmax = TPartIndex::I()->Emax();
   fNEbins = TPartIndex::I()->NEbins();
   fEGrid = TPartIndex::I()->EGrid();
   return kTRUE;
}

//___________________________________________________________________
void TEFstate::Draw(Option_t */*option*/)
{
}

