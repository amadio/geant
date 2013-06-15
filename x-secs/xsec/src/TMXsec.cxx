#include <TMXsec.h>
#include <TPartIndex.h>
#include <TMath.h>

TEXsec** TMXsec::fElements=0;

ClassImp(TMXsec)

//____________________________________________________________________________
TMXsec::TMXsec():
   fElems(0),
   fNElems(0),
   fTotXS(0), 
   fRelXS(0)
{
}

//____________________________________________________________________________
TMXsec::TMXsec(Int_t z[], Int_t /*a*/[], Float_t w[], Int_t nel, Float_t dens, Bool_t weight):
   fElems(0),
   fNElems(nel),
   fTotXS(0), 
   fRelXS(0)
{
   // Create a mixture material, we support only natural materials for the moment
   // so we ignore a (i.e. we consider it == 0)
   fElems = new Int_t[fNElems];
   if(!fElements) {
      fElements = new TEXsec*[TEXsec::NElem()];
      memset(fElements,0,TEXsec::NElem()*sizeof(TEXsec*));
   }
   for(Int_t i=0; i<fNElems; ++i)
      if((fElements[z[i]-1] = TEXsec::GetElement(z[i]))) 
	 fElems[i] = z[i]-1;
	 else Error("TMXsec","Element %d not found\n",z[i]);
   Double_t *ratios = new Double_t[fNElems];
   Double_t hnorm=0;
   for(Int_t i=0; i<fNElems; ++i) {
      ratios[i] = w[i];
      if(weight) ratios[i]/=TEXsec::WEle(z[i]);
      hnorm+=ratios[i]*TEXsec::WEle(z[i]);
   }

   if(weight) printf("By weight: ");
   else       printf("By number: ");

   for(Int_t i=0; i<fNElems; ++i) {
      ratios[i]*=TMath::Na()*1e-24*dens/hnorm;
      printf("%d %f ",z[i],ratios[i]);
   }
   printf("\n");

   // Build table with total x-sections for all mate / parts

   Int_t nbins = fElements[0]->NEbins();
   Double_t emin = fElements[0]->Emin();
   Double_t emax = fElements[0]->Emax();
   Double_t edelta = TMath::Exp(TMath::Log(emax/emin)/(nbins-1));
   Int_t totindex = TPartIndex::I()->ProcIndex("Total");
   Int_t npart = TPartIndex::I()->NPartReac();
   // Layout part1 { en<1> { tot<1>, ... , tot<fNElems>}, .....en<nbins> {tot<1>, ..., tot<fNElems>}}
   
   fRelXS = new Float_t[npart*nbins*fNElems];
   fTotXS = new Float_t[npart*nbins];
   memset(fTotXS,0,npart*nbins*sizeof(Float_t));

   for(Int_t ip=0; ip<npart; ++ip) {
      Int_t ibase = ip*(nbins*fNElems);
      Double_t en = emin;
      for(Int_t ie=0; ie<nbins; ++ie) {
	 Int_t ibin = ibase + ie*fNElems;
	 for(Int_t iel=0; iel<fNElems; ++iel) {
	    fRelXS[ibin+iel] = fElements[iel]->XS(ip,totindex,en)*ratios[iel];
	    fTotXS[ip*nbins+ie]+=fRelXS[ibin+iel];
	 }
	 if(fTotXS[ip*nbins+ie]) {
	    fTotXS[ip*nbins+ie]=1./fTotXS[ip*nbins+ie];
	    for(Int_t iel=0; iel<fNElems; ++iel) fRelXS[ibin+iel]*=fTotXS[ip*nbins+ie];
	 }
	 en*=edelta;
      }
   }
}
