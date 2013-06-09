#include "TPartIndex.h"

ClassImp(TPartIndex)

const char* TPartIndex::fPName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inHelastic",
			   "Elastic","Capture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
			   "Killer"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
			     2012,2013,2014,4131,7403};

const char* TPartIndex::fMatSymbol[NMAT]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
				"Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V",
				"Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
				"Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh",
				"Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
				"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho",
				"Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
				"Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
				"Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
				"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg",
				"Cn","Uut","Fl","Uup","Lv","Uus","Uuo"};
const char* TPartIndex::fMatName[NMAT]={"Hydrogen","Helium","Lithium","Beryllium","Boron",
				"Carbon","Nitrogen","Oxygen","Fluorine","Neon",
				"Sodium","Magnesium","Aluminium","Silicon","Phosphorus",
				"Sulfur","Chlorine","Argon","Potassium","Calcium",
				"Scandium","Titanium","Vanadium","Chromium","Manganese",
				"Iron","Cobalt","Nickel","Copper","Zinc","Gallium",
				"Germanium","Arsenic","Selenium","Bromine","Krypton",
				"Rubidium","Strontium","Yttrium","Zirconium","Niobium",
				"Molybdenum","Technetium","Ruthenium","Rhodium","Palladium",
				"Silver","Cadmium","Indium","Tin","Antimony","Tellurium",
				"Iodine","Xenon","Caesium","Barium","Lanthanum","Cerium",
				"Praseodymium","Neodymium","Promethium","Samarium",
				"Europium","Gadolinium","Terbium","Dysprosium","Holmium",
				"Erbium","Thulium","Ytterbium","Lutetium","Hafnium",
				"Tantalum","Tungsten","Rhenium","Osmium","Iridium",
				"Platinum","Gold","Mercury","Thallium","Lead","Bismuth",
				"Polonium","Astatine","Radon","Francium","Radium","Actinium",
				"Thorium","Protactinium","Uranium","Neptunium","Plutonium",
				"Americium","Curium","Berkelium","Californium","Einsteinium",
				"Fermium","Mendelevium","Nobelium","Lawrencium","Rutherfordium",
				"Dubnium","Seaborgium","Bohrium","Hassium","Meitnerium",
				"Darmstadtium","Roentgenium","Copernicium","Ununtrium",
				"Flerovium","Ununpentium","Livermorium","Ununseptium",
				"Ununoctium"};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   fNPart(0),
   fNPart30(0),
   fPDG(0),
   fPnames(0),
   fNReac(FNPREA)
{
   /*   for(Int_t i=0; i<fNProc; ++i) {
      Int_t sl = strlen(pname[i]);
      fPName[i]=new char[sl+1];
      strncpy(fPName[i],pname[i],sl);
      fPName[i][sl]='\0';
      fPCode[i]=pcode[i];
      }*/
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

//______________________________________________________________________________
void TPartIndex::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPartIndex.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TPartIndex::Class(),this);
      fgPartIndex = this;
   } else {
      R__b.WriteClassBuffer(TPartIndex::Class(),this);
   }
}




