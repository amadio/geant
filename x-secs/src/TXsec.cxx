#include <TXsec.h>
#include <TMath.h>

//___________________________________________________________________
TXsec::TXsec():
   fMat(0),
   fEmin(0),
   fEmax(0),
   fNen(0),
   fElDelta(0),
   fXsec(0),
   fNpart(0),
   fPcode(0),
   fRcode(0),
   fNXsec(0),
   fXtemp(0)
{
}

//___________________________________________________________________
TXsec::TXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen):
   fMat(z*1000+a),
   fEmin(emin),
   fEmax(emax),
   fNen(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNen),
   fXsec(0),
   fNpart(0),
   fPcode(0),
   fRcode(0),
   fNXsec(0),
   fXtemp(0)
{
}

//___________________________________________________________________
TXsec::~TXsec() {
   delete fXsec;
   delete fPcode;
   delete fRcode;
   TXtemp *ptmp = fXtemp;
   TXtemp *ttmp = 0;
   while(ptmp) {
      ttmp = ptmp;
      ptmp = ptmp->fNext;
      delete ttmp;
   }
}
