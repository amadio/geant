#include "TPFstate.h"
#include "TFinState.h"
#include "TMath.h"
#include <TFile.h>
#include <TRandom.h>

Int_t TPFstate::fVerbose=0;

ClassImp(TPFstate)

//_________________________________________________________________________
TPFstate::TPFstate():
   fPDG(0),
   fNEbins(0),
   fNReac(0),
   fNEFstat(0),
   fNFstat(0),
   fEmin(0),
   fEmax(0),
   fEilDelta(0),
   fEGrid(TPartIndex::I()->EGrid()),
   fFstat(0)
{
   Int_t np=TPartIndex::I()->NProc();
   while(np--) fRdict[np]=fRmap[np]=-1;
}

//_________________________________________________________________________
TPFstate::TPFstate(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]):
   fPDG(pdg),
   fNEbins(TPartIndex::I()->NEbins()),
   fNReac(nreac),
   fNEFstat(nfstat),
   fNFstat(fNEbins*fNEFstat),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fEilDelta((fNEbins-1)/TMath::Log(fEmax/fEmin)),
   fEGrid(TPartIndex::I()->EGrid()),
   fFstat(new TFinState[fNFstat])
{
   Int_t np=TPartIndex::I()->NProc();
   while(np--) {
      fRdict[dict[np]]=np;
      fRmap[np]=dict[np];
   }
   // consistency
   for(Int_t i=0; i<fNReac; ++i) 
      if(fRdict[fRmap[i]] != i) 
	 Fatal("SetPartXS","Dictionary mismatch for!");

}

//_________________________________________________________________________
TPFstate::~TPFstate()
{
   delete [] fFstat;
}

//_________________________________________________________________________
Bool_t TPFstate::SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]) 
{
   fPDG = pdg;
   fNEbins = TPartIndex::I()->NEbins();
   fNReac = nreac;
   fNEFstat = nfstat;
   fNFstat = fNEbins*fNReac;
   fEmin = TPartIndex::I()->Emin();
   fEmax = TPartIndex::I()->Emax();
   fEGrid = TPartIndex::I()->EGrid();
   fEilDelta = (fNEbins-1)/TMath::Log(fEmax/fEmin);
   fFstat = new TFinState[fNFstat];

   Int_t np=TPartIndex::I()->NProc();
   while(np--) {
      fRdict[dict[np]]=np;
      fRmap[np]=dict[np];
   }
   // consistency
   for(Int_t i=0; i<fNReac; ++i) 
      if(fRdict[fRmap[i]] != i) 
	 Fatal("SetPart","Dictionary mismatch for!");
   return kTRUE;
}

//_________________________________________________________________________
Bool_t TPFstate::SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[], TFinState vecfs[])
{
  fPDG = pdg;
  fNEbins = TPartIndex::I()->NEbins();
  fNReac = nreac;
  fNEFstat = nfstat;
  fNFstat = fNEbins*fNReac;
  fEmin = TPartIndex::I()->Emin();
  fEmax = TPartIndex::I()->Emax();
  fEGrid = TPartIndex::I()->EGrid();
  fEilDelta = (fNEbins-1)/TMath::Log(fEmax/fEmin);
  fFstat = vecfs;

  for(Int_t np=0; np<nreac; ++np) {
    fRdict[dict[np]]=np;
    fRmap[np]=dict[np];
  }
  // consistency
  for(Int_t i=0; i<fNReac; ++i)
    if(fRdict[fRmap[i]] != i)
      Fatal("SetPart","Dictionary mismatch for!");
  return kTRUE;
}

//_________________________________________________________________________
Bool_t TPFstate::SetFinState(Double_t en, Int_t reac, const Float_t weight[], const Float_t kerma[], 
			     const Int_t npart[], const Float_t mom[], const Int_t pid[], const Char_t surv[])
{
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
   if(en1>en || en2<en) {
      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return kFALSE;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   if(xrat<1.e-6) ++ibin;
   else if(1-xrat>1.e-6) {
      Error("SetFinState","Energy %9.2g should be %9.2g or %9.2g",en,en1,en2);
      return kFALSE;
   }
   Int_t rnumber = fRdict[reac];
   Int_t ipoint = rnumber*fNEbins + ibin;
   fFstat[ipoint].SetFinState(fNEFstat,weight,kerma,npart,mom,pid,surv);
   return kTRUE;
}     

//______________________________________________________________________________
Bool_t TPFstate::SampleReac(Double_t en, Int_t preac, Float_t& kerma, Int_t& npart,  
			   const Int_t *pid, const Float_t *mom) const
{
   Double_t eta = gRandom->Rndm();
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
   if(en1>en || en2<en) {
      Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return kFALSE;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   if(eta>xrat) ++ibin;
   Int_t rnumber = fRdict[preac];
   Int_t ipoint = rnumber*fNEbins + ibin;
   return fFstat[ipoint].SampleReac(kerma,npart,pid,mom);
}


//______________________________________________________________________________
Bool_t TPFstate::GetReac(Double_t en, Int_t preac, Int_t finstat, Float_t &kerma,
                         Int_t& npart, const Int_t *&pid, const Float_t *mom) const
{
  en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
  en=en>fEGrid[0]?en:fEGrid[0];
  Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
  ibin = ibin<fNEbins-1?ibin:fNEbins-2;
  Double_t en1 = fEGrid[ibin];
  Double_t en2 = fEGrid[ibin+1];
  if(en1>en || en2<en) {
    Error("SetFinState","Wrong bin %d in interpolation: should be %f < %f < %f\n",
          ibin, en1, en, en2);
    return kFALSE;
  }
  if(en-en1>en-en) ++ibin;
  Int_t rnumber = fRdict[preac];
  Int_t ipoint = rnumber*fNEbins + ibin;
  return fFstat[ipoint].GetReac(finstat,kerma,npart,pid,mom);
}


//______________________________________________________________________________
void TPFstate::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPFstate.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPFstate::Class(),this);
      // add the energy grid
      if(!TPartIndex::I()->EGrid()) {
	 gFile->Get("PartIndex");
      }
      fEGrid = TPartIndex::I()->EGrid();
   } else {
      R__b.WriteClassBuffer(TPFstate::Class(),this);
   }
}

//_________________________________________________________________________
void TPFstate::Print(Option_t *) const
{
   printf("Particle=%d Number of x-secs=%d between %g and %g GeV",
	  fPDG, fNReac, fEGrid[0], fEGrid[fNEbins-1]);
}

//_________________________________________________________________________
Bool_t TPFstate::Resample() {
   if(fVerbose) printf("Resampling %s from \nemin = %8.2g emacs = %8.2g, nbins = %d to \n"
		       "emin = %8.2g emacs = %8.2g, nbins = %d\n",Name(),fEmin,fEmax,
		       fNEbins,TPartIndex::I()->Emin(),TPartIndex::I()->Emax(),TPartIndex::I()->NEbins());
   // Build the original energy grid
   Double_t edelta = TMath::Exp(1/fEilDelta);
   Double_t *oGrid = new Double_t[fNEbins];
   Double_t en = fEmin;
   for(Int_t ie=0; ie<fNEbins; ++ie) {
      oGrid[ie]=en;
      en*=edelta;
   }
   // Build new arrays
   Int_t oNEbins = fNEbins;
   fNEbins = TPartIndex::I()->NEbins();
   // Temporary bins holder
   Int_t *obins = new Int_t[fNEbins];
   fNEFstat = fNEbins*fNReac;
   fEmin = TPartIndex::I()->Emin();
   fEmax = TPartIndex::I()->Emax();
   Double_t oEilDelta = fEilDelta;
   fEilDelta = TPartIndex::I()->EilDelta();
   fEGrid = TPartIndex::I()->EGrid();
   TFinState *oFstat = fFstat;
   fFstat = new TFinState[fNEFstat];

   for(Int_t ien=0; ien<fNEbins; ++ien) {
      en = fEGrid[ien];
      en=en<oGrid[oNEbins-1]?en:oGrid[oNEbins-1]*0.999;
      en=en>oGrid[0]?en:oGrid[0];
      Int_t ibin = TMath::Log(fEGrid[ien]/oGrid[0])*oEilDelta;
      ibin = ibin<oNEbins-1?ibin:oNEbins-2;
      Double_t en1 = oGrid[ibin];
      Double_t en2 = oGrid[ibin+1];
      if(en1>en || en<en) {
	 Error("Interp","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	       ibin, en1, en, en2);
      }
      Double_t xrat = (en2-en)/(en2-en1);
      if(xrat<0.5) obins[ien]=ibin;
      else obins[ien]=ibin+1;
   }

   for(Int_t ir=0; ir<fNReac; ++ir) {
      Int_t ibase = ir*fNEbins;
      for(Int_t ien=0; ien<fNEbins; ++ien) {
	 fFstat[ibase+ien] = oFstat[obins[ien]];
      }
   }

   delete [] obins;
   delete [] oFstat;
   delete [] oGrid;
   return kTRUE;
}

