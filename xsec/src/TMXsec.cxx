#include <TMXsec.h>
#include <TPartIndex.h>
#include <TMath.h>
#include <TRandom.h>

ClassImp(TMXsec)

//____________________________________________________________________________
TMXsec::TMXsec():
   fNEbins(0),
   fEmin(0),
   fEmax(0),
   fEilDelta(0),
   fEGrid(0),
   fNElems(0),
   fElems(0),
   fTotXL(0), 
   fRelXS(0)
{
}

//____________________________________________________________________________
TMXsec::TMXsec(const Char_t *name, const Char_t *title, const Int_t z[], 
	       const Int_t /*a*/[], const Float_t w[], Int_t nel, 
	       Float_t dens, Bool_t weight):
   TNamed(name,title),
   fNEbins(TPartIndex::I()->NEbins()),
   fEmin(TPartIndex::I()->Emin()),
   fEmax(TPartIndex::I()->Emax()),
   fEilDelta(TPartIndex::I()->EilDelta()),
   fEGrid(TPartIndex::I()->EGrid()),
   fNElems(nel),
   fElems(new TEXsec*[fNElems]),
   fTotXL(0), 
   fRelXS(0),
   fDEdx(0)
{
   // Create a mixture material, we support only natural materials for the moment
   // so we ignore a (i.e. we consider it == 0)
   for(Int_t i=0; i<fNElems; ++i) 
      if(z[i]) fElems[i] = TEXsec::GetElement(z[i],0);
      else if(fNElems>1) Fatal("TMXsec","Cannot have vacuum in mixtures");

   if(!z[0]) return;
      
   Double_t *ratios = new Double_t[fNElems];
   Double_t *rdedx  = new Double_t[fNElems];
   Double_t hnorm=0;
   if(fNElems>1) {
      for(Int_t i=0; i<fNElems; ++i) {
	 ratios[i] = w[i];
	 if(weight) ratios[i]/=TEXsec::WEle(z[i]);
	 hnorm+=ratios[i]*TEXsec::WEle(z[i]);
      }
   } else {
      ratios[0]=1;
      hnorm=TEXsec::WEle(z[0]);
   }

   //   if(weight) printf("By weight: ");
   //   else       printf("By number: ");

   
   for(Int_t i=0; i<fNElems; ++i) {
      rdedx[i] = ratios[i]*dens/fElems[i]->Dens();
      ratios[i]*=TMath::Na()*1e-24*dens/hnorm;
      //      printf("%d %f ",z[i],ratios[i]);
   }
   //   printf("\n");

   // Build table with total x-sections for all mate / parts

   Int_t totindex = TPartIndex::I()->ProcIndex("Total");
   Int_t npart = TPartIndex::I()->NPartReac();
   // Layout part1 { en<1> { tot<1>, ... , tot<fNElems>}, .....en<nbins> {tot<1>, ..., tot<fNElems>}}
   
   if(fNElems>1) 
      fRelXS = new Float_t[npart*fNEbins*fNElems];
   fTotXL = new Float_t[npart*fNEbins];
   memset(fTotXL,0,npart*fNEbins*sizeof(Float_t));
   fDEdx = new Float_t[npart*fNEbins];
   memset(fDEdx,0,npart*fNEbins*sizeof(Float_t));

   for(Int_t ip=0; ip<npart; ++ip) {
      Int_t ibase = ip*(fNEbins*fNElems);
      for(Int_t ie=0; ie<fNEbins; ++ie) {
	 Int_t ibin = ibase + ie*fNElems;
	 if(fNElems>1) {
	    for(Int_t iel=0; iel<fNElems; ++iel) {
	       fRelXS[ibin+iel] = fElems[iel]->XS(ip,totindex,fEGrid[ie])*ratios[iel];
	       fTotXL[ip*fNEbins+ie]+=fRelXS[ibin+iel];
	       fDEdx[ip*fNEbins+ie]+=fElems[iel]->DEdx(ip,fEGrid[ie])*rdedx[iel];
	    }
	 } else {
	    fTotXL[ip*fNEbins+ie] = fElems[0]->XS(ip,totindex,fEGrid[ie])*ratios[0];
	    fDEdx[ip*fNEbins+ie] = fElems[0]->DEdx(ip,fEGrid[ie])*rdedx[0];
	 }
	 if(fTotXL[ip*fNEbins+ie]) {
	    fTotXL[ip*fNEbins+ie]=1./fTotXL[ip*fNEbins+ie];
	    if(fNElems>1) for(Int_t iel=0; iel<fNElems; ++iel) fRelXS[ibin+iel]*=fTotXL[ip*fNEbins+ie];
	 }
      }
   }
   // cleaning up
   delete [] ratios;
   delete [] rdedx;
   delete [] fDEdx; 
}

//____________________________________________________________________________
Float_t TMXsec::Xlength(Int_t part, Float_t en) {
   if(part>=TPartIndex::I()->NPartReac() || !fTotXL) 
      return TMath::Limits<Float_t>::Max();
   else {
      en=en<=fEmax?en:fEmax;
      en=en>=fEmin?en:fEmin;
      Int_t ibin = TMath::Log(en/fEmin)*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     Double_t en1 = fEmin*TMath::Exp(fElDelta*ibin);
      //     Double_t en2 = en1*fEDelta;
      Double_t en1 = fEGrid[ibin];
      Double_t en2 = fEGrid[ibin+1];
      if(en1>en || en2<en) {
	 Error("Xlength","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	       ibin, en1, en, en2);
	 return TMath::Limits<Float_t>::Max();
      }
      Double_t xrat = (en2-en)/(en2-en1);
      return xrat*fTotXL[part*fNEbins+ibin]+(1-xrat)*fTotXL[part*fNEbins+ibin+1];
   }
}

//____________________________________________________________________________
Float_t TMXsec::DEdx(Int_t part, Float_t en) {
   if(part>=TPartIndex::I()->NPartReac() || !fDEdx) 
      return 0;
   else {
      en=en<=fEmax?en:fEmax;
      en=en>=fEmin?en:fEmin;
      Int_t ibin = TMath::Log(en/fEmin)*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     Double_t en1 = fEmin*TMath::Exp(fElDelta*ibin);
      //     Double_t en2 = en1*fEDelta;
      Double_t en1 = fEGrid[ibin];
      Double_t en2 = fEGrid[ibin+1];
      if(en1>en || en2<en) {
	 Error("DEdx","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	       ibin, en1, en, en2);
	 return TMath::Limits<Float_t>::Max();
      }
      Double_t xrat = (en2-en)/(en2-en1);
      return xrat*fDEdx[part*fNEbins+ibin]+(1-xrat)*fDEdx[part*fNEbins+ibin+1];
   }
}

//____________________________________________________________________________
TEXsec* TMXsec::SampleInt(Int_t part, Double_t en, Int_t &reac) {
   if(part>=TPartIndex::I()->NPartReac() || !fTotXL) {
      reac=-1;
      return 0;
   } else {
      en=en<=fEmax?en:fEmax;
      en=en>=fEmin?en:fEmin;
      Int_t ibin = TMath::Log(en/fEmin)*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //      Double_t en1 = fEmin*TMath::Exp(fElDelta*ibin);
      //    Double_t en2 = en1*fEDelta;
      Double_t en1 = fEGrid[ibin];
      Double_t en2 = fEGrid[ibin+1];
      if(en1>en || en2<en) {
	 Error("SampleInt","Wrong bin %d in interpolation: should be %f < %f < %f\n",
	       ibin, en1, en, en2);
	 reac=-1;
	 return 0;
      }
      Int_t iel=-1;
      if(fNElems==1) {
	 iel=0;
      } else {
	 Double_t xrat = (en2-en)/(en2-en1);
	 Double_t xnorm = 1.;
	 while(iel<0) {
	    Double_t ran = xnorm*gRandom->Rndm();
	    Double_t xsum=0;
	    for(Int_t i=0; i<fNElems; ++i) {
	       xsum+=xrat*fRelXS[ibin*fNElems+i]+(1-xrat)*fRelXS[(ibin+1)*fNElems+i];
	       if(ran<=xsum) {
		  iel = i;
		  break;
	       }
	    }
	    xnorm = xsum;
	 }
      }
      reac = fElems[iel]->SampleReac(part,en);
      return fElems[iel];
   }
}

//____________________________________________________________________________
void TMXsec::Print(Option_t*) const {
   printf("Material %s %s with %d elements\n",GetName(),GetTitle(),fNElems);
   for(Int_t i=0; i<fNElems; ++i) {
      printf("%s %s\n",fElems[i]->GetName(),fElems[i]->GetTitle());
   }
}
