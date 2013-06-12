#include <TPXsec.h>
#include <TMath.h>
#include <TRandom.h>

ClassImp(TPXsec)

//_________________________________________________________________________
TPXsec::TPXsec():
   fPDG(0),
   fNEbins(0),
   fNCbins(0),
   fNXsec(0),
   fEmin(0),
   fEmax(0),
   fEDelta(0),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0),
   fdEdx(0),
   fTotXs(0),
   fXSecs(0)
{
   memset(fRdict,0,TPartIndex::I()->NProc()*sizeof(Short_t));
   memset(fRmap,0,TPartIndex::I()->NProc()*sizeof(Short_t));
}

//_________________________________________________________________________
TPXsec::TPXsec(Int_t pdg, Int_t nen, Int_t nxsec, 
       Float_t emin, Float_t emax):
   fPDG(pdg),
   fNEbins(nen),
   fNCbins(0),
   fNXsec(nxsec),
   fTotBin(fNEbins*fNXsec),
   fEmin(emin),
   fEmax(emax),
   fEDelta(TMath::Log(fEmax/fEmin)/(fNEbins-1)),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0),
   fdEdx(0),
   fXSecs(new Float_t[fNEbins*fNXsec])
{
   memset(fRdict,0,TPartIndex::I()->NProc()*sizeof(Short_t));
   memset(fRmap,0,TPartIndex::I()->NProc()*sizeof(Short_t));
}

//_________________________________________________________________________
TPXsec::~TPXsec()
{
   delete [] fMSangle;
   delete [] fMSansig;
   delete [] fMSlength;
   delete [] fMSlensig;
   delete [] fdEdx;
   delete [] fTotXs;
   delete [] fXSecs;
}

//___________________________________________________________________
Bool_t TPXsec::SetPart(Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   fPDG = pdg;
   fNEbins = nen;
   fNXsec = nxsec;
   fTotBin = fNEbins*fNXsec;
   fEmin = emin;
   fEmax = emax;
   fEDelta = TMath::Log(fEmax/fEmin)/(fNEbins-1);
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartXS(const Float_t xsec[], const Short_t dict[]) {
   delete [] fXSecs;
   fXSecs = new Float_t[fNXsec*fNEbins];
   for(Int_t jxsec=0; jxsec<fNXsec; ++jxsec)
      for(Int_t jbin=0; jbin<fNEbins; ++jbin)
	 fXSecs[jbin*fNXsec+jxsec] = xsec[jxsec*fNEbins+jbin];
   for(Int_t i=0; i<TPartIndex::I()->NReac(); ++i) fRdict[i]=fRmap[i]=-1;
   for(Int_t i=0; i<fNXsec; ++i) {
      fRdict[TPartIndex::I()->ProcIndex(dict[i])]=i;
      fRmap[i]=TPartIndex::I()->ProcIndex(dict[i]);
   }
   // consistency
   for(Int_t i=0; i<fNXsec; ++i) 
      if(fRdict[fRmap[i]] != i) 
	 Fatal("SetPartXS","Dictionary mismatch for!");

   delete [] fTotXs;
   fTotXs = new Float_t[fNEbins];
   for(Int_t i=0; i<fNEbins; ++i) {
      fTotXs[i]=0;
      for(Int_t j=0; j<fNXsec; ++j) fTotXs[i]+=fXSecs[i*fNXsec+j];
      for(Int_t j=0; j<fNXsec; ++j) fXSecs[i*fNXsec+j]/=fTotXs[i];
   }
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartIon(const Float_t dedx[]) {
   delete [] fdEdx;
   fNCbins = fNEbins;
   fdEdx = new Float_t[fNCbins];
   memcpy(fdEdx,dedx,fNCbins*sizeof(Float_t));
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartMS(const Float_t angle[], const Float_t ansig[],
			 const Float_t length[], const Float_t lensig[]) {
   fNCbins = fNEbins;

   delete [] fMSangle;
   fMSangle = new Float_t[fNCbins];
   memcpy(fMSangle,angle,fNCbins*sizeof(Float_t));

   delete [] fMSansig;
   fMSansig = new Float_t[fNCbins];
   memcpy(fMSansig,ansig,fNCbins*sizeof(Float_t));

   delete [] fMSlength;
   fMSlength = new Float_t[fNCbins];
   memcpy(fMSlength,length,fNCbins*sizeof(Float_t));

   delete [] fMSlensig;
   fMSlensig = new Float_t[fNCbins];
   memcpy(fMSlensig,lensig,fNCbins*sizeof(Float_t));

   return kTRUE;
}

//_________________________________________________________________________
void TPXsec::Print(Option_t *) const
{
   printf("Particle=%d Number of x-secs=%d between %g and %g GeV",
	  fPDG, fNXsec, fEmin, fEmax);
}

//_________________________________________________________________________
Float_t TPXsec::DEdx(Float_t en) const {
   if(!fdEdx) return 0;
   en=en<=fEmax?en:fEmax;
   en=en>=fEmin?en:fEmin;
   Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
   Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   Double_t dedx = xrat*fdEdx[ibin]+(1-xrat)*fdEdx[ibin+1];
   /*   printf("ibin %d en1 %f en %f en2 %f xs1 %f xs2 %f xrat %f xsec %f\n",
	ibin,en1,en,en2,xs1,xs2,xrat,xsec); */
   return dedx;
}

//_________________________________________________________________________
Bool_t TPXsec::MS(Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   if(!fNCbins) {
      ang = asig = len = lsig = 0;
      return kFALSE;
   }
   en=en<=fEmax?en:fEmax;
   en=en>=fEmin?en:fEmin;
   Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
   Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   ang = xrat*fMSangle[ibin]+(1-xrat)*fMSangle[ibin+1];
   asig = xrat*fMSansig[ibin]+(1-xrat)*fMSansig[ibin+1];
   len = xrat*fMSlength[ibin]+(1-xrat)*fMSlength[ibin+1];
   lsig = xrat*fMSlensig[ibin]+(1-xrat)*fMSlensig[ibin+1];
   return kTRUE;
}

//_________________________________________________________________________
Int_t TPXsec::SampleReac(Double_t en)  const {
   en=en<=fEmax?en:fEmax;
   en=en>=fEmin?en:fEmin;
   Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
   Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);
   Double_t xnorm = 1.;
   while(1) {
      Double_t ran = xnorm*gRandom->Rndm();
      Double_t xsum=0;
      for(Int_t i=0; i<fNXsec; ++i) {
	 xsum+=xrat*fXSecs[ibin*fNXsec+i]+(1-xrat)*fXSecs[(ibin+1)*fNXsec+i];
	 if(ran<=xsum) return fRmap[i];
      }
      xnorm = xsum;
   }
}

//_________________________________________________________________________
Float_t TPXsec::XS(Short_t rindex, Float_t en) const {
   en=en<=fEmax?en:fEmax;
   en=en>=fEmin?en:fEmin;
   Int_t ibin = TMath::Log(en/fEmin)/fEDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   Double_t en1 = fEmin*TMath::Exp(fEDelta*ibin);
   Double_t en2 = fEmin*TMath::Exp(fEDelta*(ibin+1));
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);

   Double_t xtot = (xrat*fTotXs[ibin]+(1-xrat)*fTotXs[ibin+1]);
   Double_t xsec = 1;
   if(rindex<TPartIndex::I()->NReac()-1) {
      Int_t rnumber = fRdict[rindex];
      if(rnumber<0) {
	 Error("XS","No %s for %s\n",TPartIndex::I()->ReacName(rindex),
	       TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(fPDG)));
	 return -1;
      }
      xsec = xrat*fXSecs[ibin*fNXsec+rnumber]+
	 (1-xrat)*fXSecs[(ibin+1)*fNXsec+rnumber];
   }
   return xsec*xtot;
}

//_________________________________________________________________________
void TPXsec::Dump() const {
   printf("Particle %d NXsec %d emin %f emax %f NEbins %d ElDelta %f\n",
	  fPDG,fNXsec,fEmin,fEmax,fNEbins,fEDelta);
   printf("MSangle %p, MSlength %p, dEdx %p, TotXs %p, XSecs %p, Rdict %p\n",
	  fMSangle, fMSlength, fdEdx, fTotXs, fXSecs, fRdict);
}

