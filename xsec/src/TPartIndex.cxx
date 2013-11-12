#include <TPartIndex.h>
#include <TString.h>
#include <TDatabasePDG.h>
#include "TMath.h"

ClassImp(TPartIndex)

const char* TPartIndex::fPrName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inElastic",
			   "Elastic","RestCapture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture","Fission",
					"Killer","Total"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
					  2012,2013,2014,4131,4141,7403,999};

const char* TPartIndex::fEleSymbol[NELEM]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
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
const char* TPartIndex::fEleName[NELEM]={"Hydrogen","Helium","Lithium","Beryllium","Boron",
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

const Float_t TPartIndex::fWElem[NELEM]={1.008,4.0026,6.94,9.0122,10.81,12.011,14.007,15.999,
					18.998,20.180,22.990,24.305,26.982,28.085,30.974,32.06,
					35.45,39.948,39.098,40.078,44.956,47.867,50.942,51.996,
					54.938,55.845,58.933,58.693,63.546,65.38,69.723,72.63,
					74.922,78.96,79.904,83.798,85.468,87.62,88.906,91.224,
					92.906,95.96,97.91,101.07,102.91,106.42,107.87,112.41,
					114.82,118.71,121.76,127.60,126.90,131.29,132.91,137.33,
					138.91,140.12,140.91,144.24,144.91,150.36,151.96,157.25,
					158.93,162.50,164.93,167.26,168.93,173.05,174.97,178.49,
					180.95,183.84,186.21,190.23,192.22,195.08,196.97,200.59,
					204.38,207.2,208.98,208.98,209.99,222.02,223.02,226.03,
					227.03,232.04,231.04,238.03,237.05,244.06,243.06,247.07,
					247.07,251.08,252.08,257.10,258.10,259.10,262.11,265.12,
					268.13,271.13,270,277.15,276.15,281.16,280.16,285.17,
					284.18,289.19,288.19,293,294,294};

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




