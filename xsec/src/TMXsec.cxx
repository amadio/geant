#include "TMXsec.h"
#include "TPartIndex.h"
#include "TMath.h"
#include "TRandom.h"
#include "GeantTrack.h"
#include "GeantPropagator.h"
#include "GeantThreadData.h"

ClassImp(TMXsec)

//____________________________________________________________________________
TMXsec::TMXsec():
   fNEbins(0),
   fNTotXL(0),
   fNCharge(0),
   fNRelXS(0),
   fEilDelta(0),
   fEGrid(0),
   fNElems(0),
   fElems(0),
   fTotXL(0), 
   fRelXS(0),
   fDEdx(0),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0)
{
}

//____________________________________________________________________________
TMXsec::~TMXsec() {
   delete [] fElems;
   delete [] fTotXL; 
   delete [] fRelXS;
   delete [] fMSangle;
   delete [] fMSansig;
   delete [] fMSlength;
   delete [] fMSlensig;
   delete [] fDEdx;
}

//____________________________________________________________________________
TMXsec::TMXsec(const Char_t *name, const Char_t *title, const Int_t z[], 
	       const Int_t /*a*/[], const Float_t w[], Int_t nel, 
	       Float_t dens, Bool_t weight):
   TNamed(name,title),
   fNEbins(0),
   fNTotXL(0),
   fNCharge(0),
   fNRelXS(0),
   fEilDelta(TPartIndex::I()->EilDelta()),
   fEGrid(TPartIndex::I()->EGrid()),
   fNElems(0),
   fElems(0),
   fTotXL(0), 
   fRelXS(0),
   fDEdx(0),
   fMSangle(0),
   fMSansig(0),
   fMSlength(0),
   fMSlensig(0)
{
   // Create a mixture material, we support only natural materials for the moment
   // so we ignore a (i.e. we consider it == 0)

   fNElems=nel;
   fElems = new TEXsec*[fNElems];
   memset(fElems,0,fNElems*sizeof(TEXsec*));

   for(Int_t i=0; i<fNElems; ++i) 
      if(z[i]) fElems[i] = TEXsec::GetElement(z[i],0);
      else if(fNElems>1) Fatal("TMXsec","Cannot have vacuum in mixtures");

   if(!z[0]) return;

   fNEbins = TPartIndex::I()->NEbins();
   Double_t *ratios = new Double_t[fNElems];
   Double_t *rdedx  = new Double_t[fNElems];
   Double_t hnorm=0;
   if(fNElems>1) {
      for(Int_t i=0; i<fNElems; ++i) {
	 ratios[i] = w[i];
	 if(weight) ratios[i]/=TPartIndex::I()->WEle(z[i]);
	 hnorm+=ratios[i]*TPartIndex::I()->WEle(z[i]);
      }
   } else {
      ratios[0]=1;
      hnorm=TPartIndex::I()->WEle(z[0]);
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
   Int_t ncharge = TPartIndex::I()->NPartCharge();
   // Layout part1 { en<1> { tot<1>, ... , tot<fNElems>}, .....en<nbins> {tot<1>, ..., tot<fNElems>}}
   
   if(fNElems>1) {
      fNRelXS = npart*fNEbins*fNElems;
      fRelXS = new Float_t[fNRelXS];
   }
   fNTotXL = npart*fNEbins;
   fTotXL = new Float_t[fNTotXL];
   memset(fTotXL,0,fNTotXL*sizeof(Float_t));
   fNCharge = ncharge*fNEbins;
   fMSangle = new Float_t[fNCharge];
   fMSansig = new Float_t[fNCharge];
   fMSlength = new Float_t[fNCharge];
   fMSlensig = new Float_t[fNCharge];
   fDEdx = new Float_t[fNCharge];
   memset(fMSangle,0,fNCharge*sizeof(Float_t));
   memset(fMSansig,0,fNCharge*sizeof(Float_t));
   memset(fMSlength,0,fNCharge*sizeof(Float_t));
   memset(fMSlensig,0,fNCharge*sizeof(Float_t));
   memset(fDEdx,0,fNCharge*sizeof(Float_t));

   for(Int_t ip=0; ip<npart; ++ip) {
      Int_t ibase = ip*(fNEbins*fNElems);
      for(Int_t ie=0; ie<fNEbins; ++ie) {
	 Int_t ibin = ibase + ie*fNElems;
	 if(fNElems>1) {
	    for(Int_t iel=0; iel<fNElems; ++iel) {
	       fRelXS[ibin+iel] = fElems[iel]->XS(ip,totindex,fEGrid[ie])*ratios[iel];
	       fTotXL[ip*fNEbins+ie]+=fRelXS[ibin+iel];
	    }
	 } else {
	    fTotXL[ip*fNEbins+ie] = fElems[0]->XS(ip,totindex,fEGrid[ie])*ratios[0];
	 }
	 if(fTotXL[ip*fNEbins+ie]) {
	    fTotXL[ip*fNEbins+ie]=1./fTotXL[ip*fNEbins+ie];
	    if(fNElems>1) for(Int_t iel=0; iel<fNElems; ++iel) fRelXS[ibin+iel]*=fTotXL[ip*fNEbins+ie];
	 }
      }
   }
   for(Int_t ip=0; ip<ncharge; ++ip) {
      Float_t ang;
      Float_t asig;
      Float_t len;
      Float_t lsig;
      for(Int_t ie=0; ie<fNEbins; ++ie) {
	 for(Int_t iel=0; iel<fNElems; ++iel) {
	    fElems[iel]->MS(ip,fEGrid[ie],ang,asig,len,lsig);
	    fMSangle[ip*fNEbins+ie]+=ang*rdedx[iel];
	    fMSansig[ip*fNEbins+ie]+=asig*rdedx[iel];
	    fMSlength[ip*fNEbins+ie]+=len*rdedx[iel];
	    fMSlensig[ip*fNEbins+ie]+=lsig*rdedx[iel];
	    fDEdx[ip*fNEbins+ie]+=fElems[iel]->DEdx(ip,fEGrid[ie])*rdedx[iel];
	 }
      }
   }
   // cleaning up
   delete [] ratios;
   delete [] rdedx;
}


//____________________________________________________________________________
Bool_t TMXsec::Prune() {
   // Prune elements
   for(Int_t iel=0; iel<TEXsec::NLdElems(); ++iel) TEXsec::Element(iel)->Prune();
   return kTRUE;
}

//____________________________________________________________________________
Float_t TMXsec::Xlength(Int_t part, Float_t en) {
  if(part>=TPartIndex::I()->NPartReac() || !fTotXL)
    return TMath::Limits<Float_t>::Max();
  else {
    en=en<=fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
    en=en>=fEGrid[0]?en:fEGrid[0];
    Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
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
Bool_t TMXsec::Xlength_v(Int_t npart, const Int_t part[], const Float_t en[], Double_t lam[])
{
  Double_t ene;
  for(Int_t ip=0; ip<npart; ++ip) {
    if(part[ip]>=TPartIndex::I()->NPartReac() || !fTotXL)
      lam[ip]=TMath::Limits<Float_t>::Max();
    else {
      ene=en[ip]<=fEGrid[fNEbins-1]?en[ip]:fEGrid[fNEbins-1]*0.999;
      ene=en[ip]>=fEGrid[0]?en[ip]:fEGrid[0];
      Int_t ibin = TMath::Log(ene/fEGrid[0])*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     Double_t en1 = fEmin*TMath::Exp(fElDelta*ibin);
      //     Double_t en2 = en1*fEDelta;
      Double_t en1 = fEGrid[ibin];
      Double_t en2 = fEGrid[ibin+1];
      if(en1>ene || en2<ene) {
        Error("Xlength","Wrong bin %d in interpolation: should be %f < %f < %f\n",
              ibin, en1, ene, en2);
        lam[ip] = TMath::Limits<Float_t>::Max();
      }
      Double_t xrat = (en2-ene)/(en2-en1);
      lam[ip] = xrat*fTotXL[part[ip]*fNEbins+ibin]+(1-xrat)*fTotXL[part[ip]*fNEbins+ibin+1];
    }
  }
  return kTRUE;
}

//______________________________________________________________________________
void TMXsec::ProposeStep(Int_t ntracks, GeantTrack_v &tracks, Int_t tid){
// Propose step for the first ntracks in the input vector of tracks and write to
// tracks.fPstepV[]
// -should be called only for particles with reaction (first fNPartReac particle 
// in TPartIndex::fPDG[]); the case if partindex>=fNPartReac is handled in the if

   Double_t  energy;	//Ekin
   Int_t     ibin; 
   Int_t     ipart;//G5 particle index i.e. index of the particle in TPartIndex::fPDG[]
 
   // tid-based rng
   Double_t *rndArray = GeantPropagator::Instance()->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);
   
   for (Int_t i=0; i<ntracks; ++i) {
      ipart  = tracks.fG5codeV[i];
      if(ipart >= TPartIndex::I()->NPartReac() || !fTotXL) {
        tracks.fPstepV[i] = -TMath::Limits<Float_t>::Max();
      } else {
        energy = tracks.fEV[i] - tracks.fMassV[i];
        energy=energy<=fEGrid[fNEbins-1]?energy:fEGrid[fNEbins-1]*0.999;
        energy=energy>=fEGrid[0]?energy:fEGrid[0];

        ibin = TMath::Log(energy/fEGrid[0])*fEilDelta;
        ibin = ibin<fNEbins-1?ibin:fNEbins-2;

        Double_t en1 = fEGrid[ibin];
        Double_t en2 = fEGrid[ibin+1];

        Double_t xrat = (en2-energy)/(en2-en1);	
        //get interpolated -(total mean free path)
        tracks.fPstepV[i] = -(xrat*fTotXL[ipart*fNEbins+ibin]+(1-xrat)*fTotXL[ipart*fNEbins+ibin+1]);
     }
   }

   //compute proposed step based on -(total mean free path)
   for (Int_t i=0; i<ntracks; ++i) {
      tracks.fPstepV[i] *= TMath::Log(rndArray[i]);
   }

}

//get MS angles   
//____________________________________________________________________________
Float_t TMXsec::MS(Int_t ipart, Float_t energy) {
   Int_t     ibin; 
   
   if(ipart >= TPartIndex::I()->NPartCharge())
     return  0.;

   energy=energy<=fEGrid[fNEbins-1]?energy:fEGrid[fNEbins-1]*0.999;
   energy=energy>=fEGrid[0]?energy:fEGrid[0];

   ibin = TMath::Log(energy/fEGrid[0])*fEilDelta;
   ibin = ibin<fNEbins-1?ibin:fNEbins-2;

   Double_t en1 = fEGrid[ibin];
   Double_t en2 = fEGrid[ibin+1];

   if(en1>energy || en2<energy) {
     Error("MS","Wrong bin %d in interpolation: should be %f < %f < %f\n",
           ibin, en1, energy, en2);
      return  0.;
   }

   Double_t xrat = (en2-energy)/(en2-en1);
   return xrat*fMSangle[ipart*fNEbins+ibin]+(1-xrat)*fMSangle[ipart*fNEbins+ibin+1];
}

//____________________________________________________________________________
Float_t TMXsec::DEdx(Int_t part, Float_t en) {
  if(part>=TPartIndex::I()->NPartCharge() || !fDEdx)
    return 0;
  else {
    en=en<=fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
    en=en>=fEGrid[0]?en:fEGrid[0];
    Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
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
Bool_t TMXsec::DEdx_v(Int_t npart, const Int_t part[], const Float_t en[], Float_t de[]) {
  Double_t ene;
  if(!fDEdx) {
    memset(de,0,npart*sizeof(Float_t));
    return kFALSE;
  }
  for(Int_t ip=0; ip<npart; ++ip) {
    if(part[ip]>=TPartIndex::I()->NPartCharge())
      de[ip]=0;
    else {
      ene=en[ip]<=fEGrid[fNEbins-1]?en[ip]:fEGrid[fNEbins-1]*0.999;
      ene=en[ip]>=fEGrid[0]?en[ip]:fEGrid[0];
      Int_t ibin = TMath::Log(ene/fEGrid[0])*fEilDelta;
      ibin = ibin<fNEbins-1?ibin:fNEbins-2;
      //     Double_t en1 = fEmin*TMath::Exp(fElDelta*ibin);
      //     Double_t en2 = en1*fEDelta;
      Double_t en1 = fEGrid[ibin];
      Double_t en2 = fEGrid[ibin+1];
      if(en1>ene || en2<ene) {
        Error("DEdx","Wrong bin %d in interpolation: should be %f < %f < %f\n",
              ibin, en1, ene, en2);
        de[ip]=TMath::Limits<Float_t>::Max();
      }
      Double_t xrat = (en2-ene)/(en2-en1);
      de[ip] = xrat*fDEdx[part[ip]*fNEbins+ibin]+(1-xrat)*fDEdx[part[ip]*fNEbins+ibin+1];
    }
  }
  return kTRUE;
}


//____________________________________________________________________________
void TMXsec::Eloss(Int_t ntracks, GeantTrack_v &tracks)
{
// -should be called only charged particles (first fNPartCharge particle 
// in TPartIndex::fPDG[]); the case ipartindex>=fNPartCharge is handled now in the if
// Compute energy loss for the first ntracks in the input vector and update
// tracks.fEV, tracks.fPV and tracks.EdepV
   Double_t energy, dedx, edepo;
   Int_t ibin, ipart;

   Double_t energyLimit = GeantPropagator::Instance()->fEmin;
   for (Int_t i=0; i<ntracks; ++i) {
      ipart = tracks.fG5codeV[i];

      if(ipart>=TPartIndex::I()->NPartCharge())
       continue;

      energy = tracks.fEV[i] - tracks.fMassV[i]; // tabulated on kinetic
      if (energy<=fEGrid[0]) dedx = fDEdx[ipart*fNEbins];
      else if (energy>=fEGrid[fNEbins-1]) dedx = fDEdx[ipart*fNEbins+fNEbins-1]; 
      else {
         ibin = TMath::Log(energy/fEGrid[0])*fEilDelta;
         ibin = ibin<fNEbins-1?ibin:fNEbins-2;
         Double_t en1 = fEGrid[ibin];
         Double_t en2 = fEGrid[ibin+1];
         Double_t xrat = (en2-energy)/(en2-en1);
         dedx = xrat*fDEdx[ipart*fNEbins+ibin]+(1-xrat)*fDEdx[ipart*fNEbins+ibin+1];
      }
      // Update energy and momentum
      Double_t gammaold = tracks.Gamma(i);
      Double_t bgold = TMath::Sqrt((gammaold-1)*(gammaold+1));
      edepo = tracks.fStepV[i]*dedx;
      if (energy-edepo < energyLimit) {
        // Particle energy below threshold
        tracks.fEdepV[i] += energy;
        tracks.fEV[i] = tracks.fMassV[i];
        tracks.fPV[i] = 0;
        tracks.fStatusV[i] = kKilled;
        //tracks.fProcessV[i] = TPartIndex::I()->ProcIndex("kRestCapture");
        tracks.fProcessV[i] = 6; //kRestCapture will need to be called
      } else { 
        tracks.fEdepV[i] += edepo; 
        tracks.fEV[i] -= edepo;
        Double_t gammanew = tracks.Gamma(i);
        Double_t bgnew = TMath::Sqrt((gammanew-1)*(gammanew+1));
        Double_t pnorm = bgnew/bgold;
        tracks.fPV[i] *= pnorm;
      }  
   }
}   

//____________________________________________________________________________
TEXsec* TMXsec::SampleInt(Int_t part, Double_t en, Int_t &reac) {
   if(part>=TPartIndex::I()->NPartReac() || !fTotXL) {
      reac=-1;
      return 0;
   } else {
      en=en<=fEGrid[fNEbins-1]?en:fEGrid[fNEbins-1]*0.999;
      en=en>=fEGrid[0]?en:fEGrid[0];
      Int_t ibin = TMath::Log(en/fEGrid[0])*fEilDelta;
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
void TMXsec::SampleInt(Int_t ntracks, GeantTrack_v &tracksin, Int_t tid){
// -should be called only for particles with reaction (first fNPartReac particle 
// in TPartIndex::fPDG[]); the case ipartindex>=fNPartReac is handled now in the if

   Double_t energy; //Ekin
   Int_t    ibin; 
   Int_t    ipart; //G5 particle index i.e. index of the particle in TPartIndex::fPDG[]

   // tid-based rng
   Double_t *rndArray = GeantPropagator::Instance()->fThreadData[tid]->fDblArray;
   GeantPropagator::Instance()->fThreadData[tid]->fRndm->RndmArray(ntracks, rndArray);

   for (Int_t t=0; t<ntracks; ++t) {
      ipart  = tracksin.fG5codeV[t];
      if(ipart>=TPartIndex::I()->NPartReac() || !fTotXL) {
        //reac=-1;
        tracksin.fProcessV[t] = -1;//this partilce doesn't have reaction
        tracksin.fEindexV[t]  = -1;//no process -> no element selection
        continue;
      } else {
        energy = tracksin.fEV[t] - tracksin.fMassV[t];
        energy=energy<=fEGrid[fNEbins-1]?energy:fEGrid[fNEbins-1]*0.999;
        energy=energy>=fEGrid[0]?energy:fEGrid[0];

        ibin = TMath::Log(energy/fEGrid[0])*fEilDelta;
        ibin = ibin<fNEbins-1?ibin:fNEbins-2;

        Double_t en1 = fEGrid[ibin];
        Double_t en2 = fEGrid[ibin+1];
//1.
        // this material/mixture (TMXsec object) is composed of fNElems elements	 	
        Int_t iel = -1; // index of elements in TEXsec** fElems  ([fNElems]) of this 		 
        if(fNElems==1) {
	  iel=0;
        } else { //sampling one element based on the elemntal realtive X-secs that
	       //have been normalized in CTR at the Ebins!; energy interp. is 
	       //included now-> while loop is to avoid problems from interp. 
	  Double_t xrat = (en2-energy)/(en2-en1);
	  Double_t xnorm = 1.;
          while(iel<0) {
	    Double_t ran = xnorm*gRandom->Rndm();
	    Double_t xsum=0;
	    for(Int_t i=0; i<fNElems; ++i) { // simple sampling from discrete p.
	       xsum+=xrat*fRelXS[ibin*fNElems+i]+(1-xrat)*fRelXS[(ibin+1)*fNElems+i];
	       if(ran<=xsum) {
		  iel = i;
		  break;
	       }
	    }
	    xnorm = xsum;
          }
        }
        //at this point the index of the element is sampled:= iel
        //the corresponding TEXsec* is fElems[iel] 

        //sample the reaction by using the TEXsec* that corresponds to the iel-th
        //element i.e. fElems[iel] for the current particle (with particle index of
        //ipart) at the current energy; will retrun with the reaction index;
        tracksin.fProcessV[t] = fElems[iel]->SampleReac(ipart, energy, rndArray[t]); 	 
        tracksin.fEindexV[t]  = fElems[iel]->Index();//index of the selected element in TTabPhysMrg::fElemXsec[] 
      
        //INFO:
        //printf("[%d]-th partcile is %s with Ekin of %f in %s :: %s happens on %s.\n",
        //t, TPartIndex::I()->PartName(ipart), energy, this->GetName(), 
        //TPartIndex::I()->ProcName(tracksin.fProcessV[t]),fElems[iel]->GetName() );	
      }
   }
}

//____________________________________________________________________________
void TMXsec::Print(Option_t*) const {
   printf("Material %s %s with %d elements\n",GetName(),GetTitle(),fNElems);
   for(Int_t i=0; i<fNElems; ++i) {
      printf("%s %s\n",fElems[i]->GetName(),fElems[i]->GetTitle());
   }
}
