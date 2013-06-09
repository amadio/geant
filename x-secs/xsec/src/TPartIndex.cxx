#include "TPartIndex.h"

ClassImp(TPartIndex)

static const char* pname[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inHelastic",
			   "Elastic","Capture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
			   "Killer"};
static const Short_t pcode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
			     2012,2013,2014,4131,7403};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   fNProc(FNPROC),
   fNPart(0),
   fNPart30(0),
   fPDG(0),
   fPnames(0),
   fNReac(FNPREA)
{
   for(Int_t i=0; i<fNProc; ++i) {
      Int_t sl = strlen(pname[i]);
      fPName[i]=new char[sl+1];
      strncpy(fPName[i],pname[i],sl);
      fPName[i][sl]='\0';
      fPCode[i]=pcode[i];
   }
   memset(fPDGReac,0,fNReac*sizeof(Short_t));
}

//___________________________________________________________________
TPartIndex::~TPartIndex() {
}

//___________________________________________________________________
Short_t TPartIndex::ProcIndex(Short_t proccode) const {
   Short_t ip=fNProc;
   while(ip--) if(fPCode[ip]==proccode) break;
   return ip;
}

//___________________________________________________________________
const Char_t *TPartIndex::ProcNameCode(Int_t pcode) const {
   Int_t ip=ProcIndex(pcode);
   if(ip<0) return "Unknown";
   else return fPName[ip];
}

//___________________________________________________________________
const Char_t *TPartIndex::ProcNameIndex(Int_t pindex) const {
   if(pindex<0 || pindex>=fNProc) return "Unknown";
   return fPName[pindex];
}

//___________________________________________________________________
void TPartIndex::SetPartTable(char **names, Int_t *PDG, Int_t np) {
   fNPart = np;
   fNPart30 = 30*fNPart;
   fPnames = new char[fNPart30];
   fPDG = new Short_t[fNPart];
   memset(fPnames,0,fNPart30);
   for(Int_t i=0; i<fNPart; ++i) {
      fPDG[i]=PDG[i];
      strncpy(&fPnames[30*i],names[i],29);
   }
}




