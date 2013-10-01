#include "TPXsec.h"
#include "TMath.h"
#include "TRandom.h"
#include "TFile.h"

Int_t TPXsec::fVerbose=0;

ClassImp(TPXsec)

//_________________________________________________________________________
TPXsec::TPXsec():
   fPDG(0),
   fNEbins(0),
   fNCbins(0),
   fNXsec(0),
   fNTotXs(0),
   fNXSecs(0),
   fEmin(0),
   fEmax(0),
   fEilDelta(0),
   fEGrid(TPartIndex::I()->EGrid()),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0),
   fdEdx(0),
   fTotXs(0),
   fXSecs(0)
{
   Int_t np=TPartIndex::I()->NProc();
   while(np--) fRdict[np]=fRmap[np]=-1;
}

//_________________________________________________________________________
TPXsec::TPXsec(Int_t pdg, Int_t nxsec):
   fPDG(pdg),
   fNEbins(TPartIndex::I()->NEbins()),
   fNCbins(0),
   fNXsec(nxsec),
   fNTotXs(fNEbins),
   fNXSecs(fNEbins*fNXsec),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fEilDelta((fNEbins-1)/TMath::Log(fEmax/fEmin)),
   fEGrid(TPartIndex::I()->EGrid()),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0),
   fdEdx(0),
   fTotXs(0),
   fXSecs(0)
{
   Int_t np=TPartIndex::I()->NProc();
   while(np--) fRdict[np]=fRmap[np]=-1;
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

//______________________________________________________________________________
void TPXsec::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPXsec.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPXsec::Class(),this);
      // add the energy grid
      if(!TPartIndex::I()->EGrid()) {
	 gFile->Get("PartIndex");
      }
      fEGrid = TPartIndex::I()->EGrid();
   } else {
      R__b.WriteClassBuffer(TPXsec::Class(),this);
   }
}

//_________________________________________________________________________
void TPXsec::Interp(Double_t egrid[], Float_t value[], Int_t nbins, 
			   Double_t eildelta, Int_t stride, Double_t en, Float_t result[]) {
   en=en<egrid[nbins-1]?en:egrid[nbins-1]*0.999;
   en=en>egrid[0]?en:egrid[0];
   Int_t ibin = TMath::Log(en/egrid[0])*eildelta;
   ibin = ibin<nbins-1?ibin:nbins-2;
   Double_t en1 = egrid[ibin];
   Double_t en2 = egrid[ibin+1];
   if(en1>en || en2<en) {
      Error("Interp","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
   }
   Double_t xrat = (en2-en)/(en2-en1);
   for(Int_t ival=0; ival<stride; ++ival) {
      result[ival] = xrat*value[ibin*stride+ival]+
	 (1-xrat)*value[(ibin+1)*stride+ival];
   }
}

//___________________________________________________________________
Bool_t TPXsec::Prune()
{
   // 
   // Dangerous business... delete unwanted cross-sections
   //
   delete [] fTotXs;
   fTotXs=0;
   fNTotXs = 0;
   delete [] fMSangle;
   fMSangle=0;
   delete [] fMSansig;
   fMSansig=0;
   delete [] fMSlength;
   fMSlength=0;
   delete [] fMSlensig;
   fMSlensig=0;
   delete [] fdEdx;
   fdEdx=0;
   fNCbins = 0;
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::Resample() {
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
   fNTotXs = fNEbins = TPartIndex::I()->NEbins();
   if(fNCbins) fNCbins = fNEbins;
   fNXSecs = fNEbins*fNXsec;
   fEmin = TPartIndex::I()->Emin();
   fEmax = TPartIndex::I()->Emax();
   Double_t oEilDelta = fEilDelta;
   fEilDelta = TPartIndex::I()->EilDelta();
   fEGrid = TPartIndex::I()->EGrid();
   Float_t *lXSecs = new Float_t[fNXSecs];
   Float_t *lTotXs = new Float_t[fNEbins];
   // Total x-secs and partial channels
   for(Int_t ien=0; ien<fNEbins; ++ien) {
      Interp(oGrid,fXSecs,oNEbins,oEilDelta,fNXsec,fEGrid[ien],&lXSecs[ien*fNXsec]);
      Double_t xnorm=0;
      // recheck normalisation
      for(Int_t ixs=0; ixs<fNXsec; ++ixs) xnorm+=lXSecs[ien*fNXsec+ixs];
      xnorm = 1/xnorm;
      for(Int_t ixs=0; ixs<fNXsec; ++ixs) lXSecs[ien*fNXsec+ixs]*=xnorm;
      Interp(oGrid,fTotXs,oNEbins,oEilDelta,1,fEGrid[ien],&lTotXs[ien]);
   }
   delete [] fXSecs;
   fXSecs = lXSecs;
   delete [] fTotXs;
   fTotXs = lTotXs;
   // Only for charged particles
   if(fNCbins) {
      Float_t *lMSangle = new Float_t[fNCbins];
      Float_t *lMSansig = new Float_t[fNCbins];
      Float_t *lMSlength = new Float_t[fNCbins];
      Float_t *lMSlensig = new Float_t[fNCbins];
      Float_t *ldEdx = new Float_t[fNCbins];
      for(Int_t ien=0; ien<fNEbins; ++ien) {
	 Interp(oGrid,fMSangle,oNEbins,oEilDelta,1,fEGrid[ien],&lMSangle[ien]);
	 Interp(oGrid,fMSansig,oNEbins,oEilDelta,1,fEGrid[ien],&lMSansig[ien]);
	 Interp(oGrid,fMSlength,oNEbins,oEilDelta,1,fEGrid[ien],&lMSlength[ien]);
	 Interp(oGrid,fMSlensig,oNEbins,oEilDelta,1,fEGrid[ien],&lMSlensig[ien]);
	 Interp(oGrid,fdEdx,oNEbins,oEilDelta,1,fEGrid[ien],&ldEdx[ien]);
      }
      delete [] fMSangle;
      fMSangle = lMSangle;
      delete [] fMSansig;
      fMSansig = lMSansig;
      delete [] fMSlength;
      fMSlength = lMSlength;
      delete [] fMSlensig;
      fMSlensig = lMSlensig;
      delete [] fdEdx;
      fdEdx = ldEdx;
   }
   delete [] oGrid;
   return kTRUE;
}


//___________________________________________________________________
Bool_t TPXsec::SetPart(Int_t pdg, Int_t nxsec) {
   fPDG = pdg;
   fNEbins = TPartIndex::I()->NEbins();
   fNXsec = nxsec;
   fNTotXs = fNEbins;
   fNXSecs = fNEbins*fNXsec;
   fEmin = TPartIndex::I()->Emin();
   fEmax = TPartIndex::I()->Emax();
   fEGrid = TPartIndex::I()->EGrid();
   fEilDelta = TPartIndex::I()->EilDelta();
   return kTRUE;
}

//___________________________________________________________________
Bool_t TPXsec::SetPartXS(const Float_t xsec[], const Int_t dict[]) {
   delete [] fXSecs;
   fXSecs = new Float_t[fNXSecs];
   for(Int_t jxsec=0; jxsec<fNXsec; ++jxsec)
      for(Int_t jbin=0; jbin<fNEbins; ++jbin)
	 fXSecs[jbin*fNXsec+jxsec] = xsec[jxsec*fNEbins+jbin];
   for(Int_t i=0; i<TPartIndex::I()->NProc(); ++i) fRdict[i]=fRmap[i]=-1;
   for(Int_t i=0; i<fNXsec; ++i) {
      fRdict[dict[i]]=i;
      fRmap[i]=dict[i];
   }
   // consistency
   for(Int_t i=0; i<fNXsec; ++i) 
      if(fRdict[fRmap[i]] != i) 
	 Fatal("SetPartXS","Dictionary mismatch for!");

   delete [] fTotXs;
   fTotXs = new Float_t[fNTotXs];
   for(Int_t i=0; i<fNEbins; ++i) {
      fTotXs[i]=0;
      for(Int_t j=0; j<fNXsec; ++j) fTotXs[i]+=fXSecs[i*fNXsec+j];
      if(fTotXs[i]) 
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
	  fPDG, fNXsec, fEGrid[0], fEGrid[fNEbins-1]);
}

//_________________________________________________________________________
Float_t TPXsec::DEdx(Double_t en) const {
   if(!fdEdx) return 0;
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   //   Double_t en1 = fEmin*TMath::Exp(ibin/fEilDelta);
   //   Double_t en2 = fEmin*TMath::Exp((ibin+1)/fEilDelta);
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
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
Bool_t TPXsec::MS(Double_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const {
   if(!fNCbins) {
      ang = asig = len = lsig = 0;
      return kFALSE;
   }
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   //   Double_t en1 = fEmin*TMath::Exp(ibin/fEilDelta);
   //  Double_t en2 = fEmin*TMath::Exp((ibin+1)/fEilDelta);
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
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
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
   //   Double_t en1 = fEmin*TMath::Exp(ibin/fEilDelta);
   //Double_t en2 = fEmin*TMath::Exp((ibin+1)/fEilDelta);
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
   if(en1>en || en2<en) {
      Error("SampleReac","Wrong bin %d in interpolation: should be %f < %f < %f\n",
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
Bool_t TPXsec::XS_v(Int_t npart, Int_t rindex, const Double_t en[], Double_t lam[]) const
{
  //   printf("fEGrid %p\n",fEGrid);
  Double_t ene;
  for(Int_t ip=0; ip<npart; ++ip) {
    ene=en[ip]<fEGrid[fNEbins-1]?en[ip]:fEGrid[fNEbins-1]*0.999;
    ene=en[ip]>fEGrid[0]?en[ip]:fEGrid[0];
    Int_t ibin = TMath::Log(ene/fEGrid[0])*fEilDelta;
    ibin = ibin<fNEbins-1?ibin:fNEbins-2;
    //   Double_t en1 = fEmin*TMath::Exp(ibin/fEilDelta);
    //   Double_t en2 = fEmin*TMath::Exp((ibin+1)/fEilDelta);
    Double_t en1 = fEGrid[ibin];
    Double_t en2 = fEGrid[ibin+1];
    if(en1>ene || en2<ene) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
            ibin, en1, ene, en2);
      return 0;
    }
    Double_t xrat = (en2-ene)/(en2-en1);
    
    Double_t xtot = (xrat*fTotXs[ibin]+(1-xrat)*fTotXs[ibin+1]);
    Double_t xsec = 1;
    if(rindex<TPartIndex::I()->NProc()-1) {
      Int_t rnumber = fRdict[rindex];
      if(rnumber<0) {
        Error("XS","No %s for %s\n",TPartIndex::I()->ProcName(rindex),
              TPartIndex::I()->PartName(TPartIndex::I()->PartIndex(fPDG)));
        return -1;
      }
      xsec = xrat*fXSecs[ibin*fNXsec+rnumber]+
      (1-xrat)*fXSecs[(ibin+1)*fNXsec+rnumber];
    }
    lam[ip] = xsec*xtot;
  }
  return kTRUE;
}

//_________________________________________________________________________
Float_t TPXsec::XS(Int_t rindex, Double_t en) const {
   //   printf("fEGrid %p\n",fEGrid);
   en=en<fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
   en=en>fEGrid[0]?en:fEGrid[0];
   Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;
//   Double_t en1 = fEmin*TMath::Exp(ibin/fEilDelta);
//   Double_t en2 = fEmin*TMath::Exp((ibin+1)/fEilDelta);
   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];
   if(en1>en || en2<en) {
      Error("XS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	    ibin, en1, en, en2);
      return 0;
   }
   Double_t xrat = (en2-en)/(en2-en1);

   Double_t xtot = (xrat*fTotXs[ibin]+(1-xrat)*fTotXs[ibin+1]);
   Double_t xsec = 1;
   if(rindex<TPartIndex::I()->NProc()-1) {
      Int_t rnumber = fRdict[rindex];
      if(rnumber<0) {
	 Error("XS","No %s for %s\n",TPartIndex::I()->ProcName(rindex),
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
	  fPDG,fNXsec,fEGrid[0],fEGrid[fNEbins-1],fNEbins,fEilDelta);
   printf("MSangle %p, MSlength %p, dEdx %p, TotXs %p, XSecs %p, Rdict %p\n",
	  (void*)fMSangle, (void*)fMSlength, (void*)fdEdx, (void*)fTotXs, 
	  (void*)fXSecs, (const void*)fRdict);
}

