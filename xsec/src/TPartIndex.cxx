#include <TPartIndex.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include <TMath.h>

ClassImp(TPartIndex)

const char* TPartIndex::fPrName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inElastic",
			   "Elastic","RestCapture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
					"Killer","Total"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
					  2012,2013,2014,4131,7403,999};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   TNamed("PartIndex","Translation tables"),
   fNPart(0),
   fPDG(0),
   fNpReac(0),
   fNpCharge(0),
   fNEbins(0),
   fEilDelta(0),
   fEGrid(0),
   fDBPdg(TDatabasePDG::Instance())
{ 
}

//___________________________________________________________________
TPartIndex::~TPartIndex() {
   delete [] fPDG;
   delete [] fEGrid;
   delete fDBPdg;
   fgPartIndex=0;
}


//___________________________________________________________________
void TPartIndex::SetEnergyGrid(Double_t emin, Double_t emax, Int_t nbins) {
   fNEbins = nbins;
   fEilDelta = (fNEbins-1)/TMath::Log(emax/emin);
   delete [] fEGrid;
   fEGrid = new Double_t[fNEbins];
   Double_t en=emin;
   Double_t edelta=TMath::Exp(1/fEilDelta);
   for(Int_t i=0; i<fNEbins; ++i) {
      fEGrid[i]=en;
      en*=edelta;
   }
}

//___________________________________________________________________
Int_t TPartIndex::ProcIndex(Int_t proccode) const {
   Short_t ip=fNProc;
   while(ip--) if(fPCode[ip]==proccode) break;
   return ip;
}

//___________________________________________________________________
const Char_t *TPartIndex::ProcName(Int_t proc) const {
   if(proc<0 || proc>=fNProc) return "Unknown";
   return fPrName[proc];
}

//___________________________________________________________________
void TPartIndex::SetPartTable(const Int_t *vpdg, Int_t np) {
   fNPart = np;
   delete [] fPDG;
   fPDG = new Int_t[fNPart];
   for(Int_t i=0; i<fNPart; ++i) 
      fPDG[i]=vpdg[i];
}

//______________________________________________________________________________
Int_t TPartIndex::PDG(const Char_t* pname) const {Int_t nr=fNPart;
   while(nr--) if(!strcmp(pname,fDBPdg->GetParticle(fPDG[nr])->GetName())) return fPDG[nr];
   return -12345678;
}

//______________________________________________________________________________
void TPartIndex::Print(Option_t *option) const
{
   Char_t line[120];
   TString opt = option;
   opt.ToLower();
   if(opt.Contains("particles")) {
      printf("Available particles:\n");
      memset(line,0,120);
      for(Int_t i=0; i<fNPart; ++i) {
	 const char *name = fDBPdg->GetParticle(fPDG[i])->GetName();
	 if(strlen(line)+strlen(name)+1>119) {
	    printf("%s\n",line);
	    memset(line,0,120);
	 }
	 strcat(line,name);
	 strcat(line," ");
      }
      if(strlen(line)) printf("%s\n",line);
   }
   if(opt.Contains("reactions")) {
      printf("Available reactions:\n");
      memset(line,0,120);
      strcat(line,"Total ");
      for(Int_t i=0; i<fNProc; ++i) {
	 if(strlen(line)+strlen(fPrName[i])+1>119) {
	    printf("%s\n",line);
	    memset(line,0,120);
	 }
	 strcat(line,fPrName[i]);
	 strcat(line," ");
      }
      if(strlen(line)) printf("%s\n",line);
   }
}

//______________________________________________________________________________
void TPartIndex::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPartIndex.

   if (R__b.IsReading()) {
      delete fDBPdg;
      fDBPdg = 0;
      R__b.ReadClassBuffer(TPartIndex::Class(),this);
      fgPartIndex = this;
   } else {
      R__b.WriteClassBuffer(TPartIndex::Class(),this);
   }
}




