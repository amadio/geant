#include <TPXsec.h>
#include <TMath.h>

//_________________________________________________________________________
TPXsec::TPXsec():
   fPDG(0),
   fNen(0),
   fNXsec(0),
   fEmin(0),
   fEmax(0),
   fEDelta(0),
   fMSangle(0),
   fMSlength(0),
   fdEdx(0),
   fTotXs(0),
   fXSecs(0),
   fRdict(0)
{
}

//_________________________________________________________________________
TPXsec::TPXsec(Int_t pdg, Int_t nen, Int_t nxsec, 
       Float_t emin, Float_t emax):
   fPDG(pdg),
   fNen(nen),
   fNXsec(nxsec),
   fEmin(emin),
   fEmax(emax),
   fEDelta(TMath::Log(emax/emin)/fNen),
   fMSangle(0),
   fMSlength(0),
   fdEdx(0),
   fTotXs(new Float_t[fNen]),
   fXSecs(0),
   fRdict(0)
{
}

//_________________________________________________________________________
TPXsec::~TPXsec()
{
   delete fMSangle;
   delete fMSlength;
   delete fdEdx;
   delete fTotXs;
   delete fXSecs;
   delete fRdict;
}

//_________________________________________________________________________
void TPXsec::Print(Option_t *) const
{
   printf("Particle=%d Number of x-secs=%d between %g and %g GeV",
	  fPDG, fNXsec, fEmin, fEmax);
}


