#include <TPXsec.h>
#include <TMath.h>

ClassImp(TPXsec)

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
   fTotBin(fNen*fNXsec),
   fEmin(emin),
   fEmax(emax),
   fEDelta(TMath::Log(emax/emin)/fNen),
   fMSangle(new Float_t[fNen]),
   fMSlength(new Float_t[fNen]),
   fdEdx(new Float_t[fNen]),
   fTotXs(new Float_t[fNen]),
   fXSecs(new Float_t[fNen*fNXsec]),
   fRdict(new Short_t[nxsec])
{
}

//_________________________________________________________________________
TPXsec::~TPXsec()
{
   delete [] fMSangle;
   delete [] fMSlength;
   delete [] fdEdx;
   delete [] fTotXs;
   delete [] fXSecs;
   delete [] fRdict;
}

//___________________________________________________________________
Bool_t TPXsec::SetPart(Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   fPDG = pdg;
   fNen = nen;
   fNXsec = nxsec;
   fTotBin = fNen*fNXsec;
   fEmin = emin;
   fEmax = emax;
   fEDelta = TMath::Log(emax/emin)/fNen;
   fMSangle = new Float_t[fNen];
   fMSlength = new Float_t[fNen];
   fdEdx = new Float_t[fNen];
   fTotXs = new Float_t[fNen];
   fXSecs = new Float_t[fNen*fNXsec];
   fRdict = new Short_t[nxsec];
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartXS(const Float_t xsec[], const Short_t dict[]) {
   delete [] fXSecs;
   fXSecs = new Float_t[fNXsec*fNen];
   memcpy(fXSecs,xsec,fNXsec*fNen*sizeof(Float_t));
   delete [] fRdict;
   fRdict = new Short_t[fNXsec];
   memcpy(fRdict,dict,fNXsec*sizeof(Short_t));
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartIon(const Float_t dedx[]) {
   delete [] fdEdx;
   fdEdx = new Float_t[fNen];
   memcpy(fdEdx,dedx,fNen*sizeof(Float_t));
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::Finalise() {
   delete [] fTotXs;
   fTotXs = new Float_t[fNen];
   for(Int_t i=0; i<fNen; ++i) {
      fTotXs[i]=0;
      for(Int_t j=0; j<fNXsec; ++j) fTotXs[i]+=fXSecs[j*fNen+i];
   }
   return kTRUE;
}

//_________________________________________________________________________
void TPXsec::Print(Option_t *) const
{
   printf("Particle=%d Number of x-secs=%d between %g and %g GeV",
	  fPDG, fNXsec, fEmin, fEmax);
}

//_________________________________________________________________________
Float_t TPXsec::XS(Short_t rcode, Float_t en) const {
   for(Int_t i=0; i<fNXsec; ++i) 
      if(rcode == fRdict[i]) {
	 if(en>=fEmax) return fXSecs[fNen-1];
	 if(en<=fEmin) return fXSecs[0];
	 Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
	 Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
	 Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
	 Double_t xs1 = fXSecs[i*fNen+ibin];
	 Double_t xs2 = fXSecs[i*fNen+ibin+1];
	 printf("ibin %d %f < %f < %f\n",ibin,en1,en,en2);
	 return (en2-en)*xs1/(en2-en1)+(en-en1)*xs2/(en2-en1);
      }
   Error("XS","No reaction code %d\n",rcode);
   return 0;
}

//_________________________________________________________________________
void TPXsec::Dump() const {
   printf("Particle %d NXsec %d emin %f emax %f NEbins %d ElDelta %f\n",
	  fPDG,fNXsec,fEmin,fEmax,fNen,fEDelta);
   printf("MSangle %p, MSlength %p, dEdx %p, TotXs %p, XSecs %p, Rdict %p\n",
	  fMSangle, fMSlength, fdEdx, fTotXs, fXSecs, fRdict);
}

