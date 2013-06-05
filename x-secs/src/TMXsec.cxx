#include <TMXsec.h>
#include <TPXsec.h>
#include <TMath.h>

//___________________________________________________________________
TMXsec::TMXsec():
   fMat(0),
   fEmin(0),
   fEmax(0),
   fNen(0),
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
   fNen(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNen),
   fNpart(np),
   fPXsec(new TPXsec[fNpart]),
   fCuts(0)
{
}

//___________________________________________________________________
TMXsec::~TMXsec() {
   delete fPXsec;
   delete fCuts;
}
