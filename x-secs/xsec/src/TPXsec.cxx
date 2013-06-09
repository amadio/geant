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
   fXSecs(0)
{
   memset(fRdict,0,TPartIndex::I()->NProc()*sizeof(Short_t));
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
   fEDelta(TMath::Log(fEmax/fEmin)/(fNen-1)),
   fMSangle(new Float_t[fNen]),
   fMSlength(new Float_t[fNen]),
   fdEdx(new Float_t[fNen]),
   fTotXs(new Float_t[fNen]),
   fXSecs(new Float_t[fNen*fNXsec])
{
   memset(fRdict,0,TPartIndex::I()->NProc()*sizeof(Short_t));
}

//_________________________________________________________________________
TPXsec::~TPXsec()
{
   delete [] fMSangle;
   delete [] fMSlength;
   delete [] fdEdx;
   delete [] fTotXs;
   delete [] fXSecs;
}

//___________________________________________________________________
Bool_t TPXsec::SetPart(Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   fPDG = pdg;
   fNen = nen;
   fNXsec = nxsec;
   fTotBin = fNen*fNXsec;
   fEmin = emin;
   fEmax = emax;
   fEDelta = TMath::Log(fEmax/fEmin)/(fNen-1);
   fMSangle = new Float_t[fNen];
   fMSlength = new Float_t[fNen];
   fdEdx = new Float_t[fNen];
   fTotXs = new Float_t[fNen];
   fXSecs = new Float_t[fNen*fNXsec];
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartXS(const Float_t xsec[], const Short_t dict[]) {
   delete [] fXSecs;
   fXSecs = new Float_t[fNXsec*fNen];
   memcpy(fXSecs,xsec,fNXsec*fNen*sizeof(Float_t));
   for(Int_t i=0; i<fNXsec; ++i) fRdict[0]=-1;
   memset(fRdict,0,TPartIndex::I()->NProc()*sizeof(Short_t));
   for(Int_t i=0; i<fNXsec; ++i) fRdict[i]=TPartIndex::I()->ProcIndex(dict[i]);
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
Float_t TPXsec::XS(Short_t rindex, Float_t en) const {
   en=en<=fEmax?en:fEmax;
   en=en>=fEmin?en:fEmin;
   Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
   ibin = ibin<fNen-1?ibin:fNen-2;
   Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
   Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
   Double_t xs1 = fXSecs[rindex*fNen+ibin];
   Double_t xs2 = fXSecs[rindex*fNen+ibin+1];
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   Double_t xsec = xrat*xs1+(1-xrat)*xs2;
   printf("ibin %d en1 %f en %f en2 %f xs1 %f xs2 %f xrat %f xsec %f\n",
	  ibin,en1,en,en2,xs1,xs2,xrat,xsec);
   return xsec;
}

//_________________________________________________________________________
void TPXsec::Dump() const {
   printf("Particle %d NXsec %d emin %f emax %f NEbins %d ElDelta %f\n",
	  fPDG,fNXsec,fEmin,fEmax,fNen,fEDelta);
   printf("MSangle %p, MSlength %p, dEdx %p, TotXs %p, XSecs %p, Rdict %p\n",
	  fMSangle, fMSlength, fdEdx, fTotXs, fXSecs, fRdict);
}

