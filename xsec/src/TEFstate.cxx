#include "TAxis.h"
#include "TCanvas.h"
#include "TEFstate.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
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
   fNRpart(0),
   fPFstate(0)
{
}

//___________________________________________________________________
TEFstate::TEFstate(Int_t z, Int_t a, Float_t dens):
   TNamed(TPartIndex::I()->EleSymb(z),TPartIndex::I()->EleName(z)),
   fEle(z*10000+a*10),
   fDens(dens),
   fAtcm3(fDens*TMath::Na()*1e-24/TPartIndex::I()->WEle(z)),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fNEbins(TPartIndex::I()->NEbins()),
   fEilDelta(TPartIndex::I()->EilDelta()),
   fEGrid(TPartIndex::I()->EGrid()),
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
  return fPFstate[kpart].SetPart(pdg,nfstat,nreac,dict,vecfs);
}

//___________________________________________________________________
Bool_t TEFstate::AddPartFS(Int_t kpart, Double_t en, Int_t reac, const Float_t weight[], 
			   const Float_t kerma[], const Int_t npart[], const Float_t mom[], 
			   const Int_t pid[], const Char_t surv[])
{
   return fPFstate[kpart].SetFinState(en,reac,weight,kerma,npart,mom,pid,surv);
}

//___________________________________________________________________
Bool_t TEFstate::SampleReac(Int_t pindex, Double_t en, Int_t preac,
			 Float_t& kerma, Int_t &npart, const Int_t *pid, const Float_t *mom) const {
   return fPFstate[pindex].SampleReac(preac,en,kerma,npart,pid,mom);
}

//___________________________________________________________________
Bool_t TEFstate::GetReac(Int_t pindex, Double_t en, Int_t preac, Int_t ifs,
                              Float_t& kerma, Int_t &npart, const Int_t *&pid, const Float_t *mom) const
{
  return fPFstate[pindex].GetReac(en,preac,ifs,kerma,npart,pid,mom);
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
      if(FloatDiff(TPartIndex::I()->Emin(),fElements[fNLdElems]->Emin(),1e-7) ||
	 FloatDiff(TPartIndex::I()->Emax(),fElements[fNLdElems]->Emax(),1e-7) ||
	 TPartIndex::I()->NEbins() != fElements[fNLdElems]->NEbins())
	 // we have to resize the energy grid of the element
	 fElements[fNLdElems]->Resample();
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

