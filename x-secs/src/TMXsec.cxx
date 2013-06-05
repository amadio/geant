#include <TMXsec.h>
#include <TPXsec.h>
#include <TMath.h>

ClassImp(TMXsec)

//___________________________________________________________________
TMXsec::TMXsec():
   fMat(0),
   fEmin(0),
   fEmax(0),
   fNEbins(0),
   fElDelta(0),
   fNpart(0),
   fPXsec(0),
   fCuts(0)
{
}

//___________________________________________________________________
TMXsec::TMXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np):
   fMat(z*10000+a*10),
   fEmin(emin),
   fEmax(emax),
   fNEbins(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNEbins),
   fNpart(np),
   fPXsec(new TPXsec[fNpart]),
   fCuts(new Double_t[fNpart])
{
}

//___________________________________________________________________
TMXsec::~TMXsec() {
   delete fPXsec;
   delete fCuts;
}

//___________________________________________________________________
Bool_t TMXsec::AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   return fPXsec[kpart].SetPart(pdg,nen,nxsec,emin,emax);
}


//___________________________________________________________________
Bool_t TMXsec::AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]) {
   return fPXsec[kpart].SetPartXS(xsec,dict);
}

//___________________________________________________________________
Bool_t TMXsec::AddPartIon(Int_t kpart, const Float_t dedx[]) {
   return fPXsec[kpart].SetPartIon(dedx);
}

//___________________________________________________________________
Bool_t TMXsec::Finalise() {
   Bool_t retok = kTRUE;
   for(Int_t i=0; i<fNpart; ++i) 
      retok &= fPXsec[i].Finalise();
   return retok;
}

//___________________________________________________________________
Float_t TMXsec::XS(Int_t pdg, Short_t rcode, Float_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].XS(rcode,en);
   return -1;
}

//___________________________________________________________________
void TMXsec::Dump() const {
   printf("Material %d emin %f emax %f NEbins %d ElDelta %f Npart %d\n",
	  fMat,fEmin,fEmax,fNEbins,fElDelta,fNpart);
   for(Int_t i=0; i<fNpart; ++i) 
      fPXsec[i].Dump();
}

