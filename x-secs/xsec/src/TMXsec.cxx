#include <TMXsec.h>
#include <TPXsec.h>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TObjArray.h>
#include <TPartIndex.h>
#include <TString.h>
#include <TObjString.h>
#include <TFile.h>

ClassImp(TMXsec)

//___________________________________________________________________
TMXsec::TMXsec():
   fMat(0),
   fEmin(0),
   fEmax(0),
   fNEbins(0),
   fElDelta(0),
   fNpart(0),
   fPXsec(0),
   fCuts(0)
{
}

//___________________________________________________________________
TMXsec::TMXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np):
   TNamed(TPartIndex::MatSymb(z),TPartIndex::MatName(z)),
   fMat(z*10000+a*10),
   fEmin(emin),
   fEmax(emax),
   fNEbins(nen),
   fElDelta(TMath::Log(fEmax/fEmin)/fNEbins),
   fNpart(np),
   fPXsec(new TPXsec[fNpart]),
   fCuts(new Double_t[fNpart])
{
}

//___________________________________________________________________
TMXsec::~TMXsec() {
   delete fPXsec;
   delete fCuts;
}

//___________________________________________________________________
Bool_t TMXsec::AddPart(Int_t kpart, Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax) {
   return fPXsec[kpart].SetPart(pdg,nen,nxsec,emin,emax);
}


//___________________________________________________________________
Bool_t TMXsec::AddPartXS(Int_t kpart, const Float_t xsec[], const Short_t dict[]) {
   return fPXsec[kpart].SetPartXS(xsec,dict);
}

//___________________________________________________________________
Bool_t TMXsec::AddPartIon(Int_t kpart, const Float_t dedx[]) {
   return fPXsec[kpart].SetPartIon(dedx);
}

//___________________________________________________________________
Bool_t TMXsec::Finalise() {
   Bool_t retok = kTRUE;
   for(Int_t i=0; i<fNpart; ++i) 
      retok &= fPXsec[i].Finalise();
   return retok;
}

//___________________________________________________________________
Float_t TMXsec::XS(Int_t pdg, Short_t rcode, Float_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].XS(TPartIndex::I()->ProcIndex(rcode),en);
   return -1;
}

//___________________________________________________________________
Float_t TMXsec::XSindex(Int_t pindex, Short_t rindex, Float_t en) const {
   return fPXsec[pindex].XS(rindex,en);
}

//___________________________________________________________________
Float_t TMXsec::DEdx(Int_t pdg, Float_t en) const {
   for(Int_t i=0; i<fNpart; ++i) 
      if(pdg == fPXsec[i].PDG()) 
	 return fPXsec[i].DEdx(en);
   return -1;
}

//___________________________________________________________________
Float_t TMXsec::DEdxIndex(Int_t pindex, Float_t en) const {
   return fPXsec[pindex].DEdx(en);
}

//___________________________________________________________________
void TMXsec::DumpPointers() const {
   printf("Material %d emin %f emax %f NEbins %d ElDelta %f Npart %d\n",
	  fMat,fEmin,fEmax,fNEbins,fElDelta,fNpart);
   for(Int_t i=0; i<fNpart; ++i) 
      fPXsec[i].Dump();
}

//___________________________________________________________________
TGraph* TMXsec::XSGraph(const char* part, const char *reac, 
			Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *xsec = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t rcode = TPartIndex::I()->Rcode(reac);
   if(rcode<0) {
      Error("XSGraph","Unknown reaction %s\n",reac);
      return 0;
   }
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      xsec[i] = XS(pdg,rcode,en);
      en*=delta;
   }
   TGraph *tg = new TGraph(nbin,energy,xsec);
   memset(title,0,200);
   snprintf(title,199,"%s %s on %s",part,reac,GetTitle());
   tg->SetTitle(title);
   delete [] xsec;
   delete [] energy;
   return tg;
}

//___________________________________________________________________
TGraph* TMXsec::DEdxGraph(const char* part, 
			  Float_t emin, Float_t emax, Int_t nbin) const 
{
   Char_t title[200];
   const Double_t delta = TMath::Exp(TMath::Log(emax/emin)/(nbin-1));
   Float_t *dedx = new Float_t[nbin];
   Float_t *energy = new Float_t[nbin];
   Double_t en=emin;
   Int_t pdg = TPartIndex::I()->PDG(part);
   if(pdg == -12345678) {
      Error("XSGraph","Unknown particle %s\n",part);
      return 0;
   }
   for(Int_t i=0; i<nbin; ++i) {
      energy[i] = en;
      dedx[i] = DEdx(pdg,en);
      en*=delta;
   }
   TGraph *tg = new TGraph(nbin,energy,dedx);
   memset(title,0,200);
   snprintf(title,199,"%s dEdx on %s",part,GetTitle());
   tg->SetTitle(title);
   delete [] dedx;
   delete [] energy;
   return tg;
}

//___________________________________________________________________
void TMXsec::Draw(Option_t *option)
{
   // Draw cross sections and other physics quantities for this material
   //
   // Format is "particle,reaction,emin,emax,nbin"
   // Available reactions are: 
   // Transport,MultScatt,Ionisation,Decay,inElastic,Elastic,Capture,Brehms,PairProd",Annihilation,
   // CoulombScatt,Photoel,Compton,Conversion,Capture,Killer
   Char_t title[200];
   TString opt = option;
   TObjArray *token = opt.Tokenize(",");
   Int_t narg = token->GetEntries();
   if(narg<2) {
      Error("Draw","Must speficy at least particle and reaction");
      return;
   }
   if(gFile) gFile->Get("PartIndex");
   const Char_t *part = ((TObjString *) token->At(0))->GetName();
   const Char_t *reac = ((TObjString*) token->At(1))->GetName();
   Float_t emin=1e-3;
   if(narg>2) sscanf(((TObjString*) token->At(2))->GetName(),"%f",&emin);
   Float_t emax=1e6;
   if(narg>3) sscanf(((TObjString*) token->At(3))->GetName(),"%f",&emax);
   Int_t nbin=100;
   if(narg>4) sscanf(((TObjString*) token->At(4))->GetName(),"%d",&nbin);
   if(TPartIndex::I()->PDG(part) == -12345678) {
      Error("Draw","Particle %s does not exist\n",part);
      TPartIndex::I()->Print("particles");
      return;
   }
   TGraph *tg=0;
   if(!strcmp(reac,"dEdx")) {
      tg = DEdxGraph(part, emin, emax, nbin);
      snprintf(title,199,"%s dEdx on %s",part,GetTitle());
   } else {
      if(TPartIndex::I()->Rcode(reac)<0) {
	 Error("Draw","Reaction %s does not exist\n",reac);
	 TPartIndex::I()->Print("reactions");
	 printf("dEdx\n");
	 return;
      }
      tg = XSGraph(part, reac, emin, emax, nbin);
      snprintf(title,199,"%s %s on %s",part,reac,GetTitle());
   }
   TCanvas *tc = new TCanvas(part,title,600,400);
   tc->SetLogx();
   tg->Draw("ACL");
}

