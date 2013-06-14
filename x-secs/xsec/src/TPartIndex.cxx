#include "TPartIndex.h"
#include "TString.h"

ClassImp(TPartIndex)

const char* TPartIndex::fPrName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inElastic",
			   "Elastic","Capture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
					"Killer","Total"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
					  2012,2013,2014,4131,7403,999};

const char* TPartIndex::fEleSymbol[NMAT]={"H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
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
const char* TPartIndex::fEleName[NMAT]={"Hydrogen","Helium","Lithium","Beryllium","Boron",
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

const Float_t TPartIndex::fWElem[NMAT]={1.008,4.0026,6.94,9.0122,10.81,12.011,14.007,15.999,
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


static const Int_t kPDG[FNPART]={
-1000020040,-1000020030,-1000010030,-1000010020,-100012210,-100012110,-100002210,-100002110,-9000211,-100325,
-100323,-100321,-100315,-100313,-100311,-100213,-100211,-53122,-52214,-52114,-43122,-42212,-42124,-42112,-41214,
-33324,-33314,-33122,-32224,-32214,-32212,-32124,-32114,-32112,-31214,-31114,-30323,-30313,-30213,-23324,-23314,
-23224,-23222,-23214,-23212,-23126,-23124,-23122,-23114,-23112,-22224,-22222,-22214,-22212,-22124,-22122,-22114,
-22112,-21214,-21212,-21114,-21112,-20323,-20313,-20213,-13326,-13324,-13316,-13314,-13226,-13224,-13222,-13216,
-13214,-13212,-13126,-13124,-13122,-13116,-13114,-13112,-12226,-12224,-12222,-12218,-12216,-12214,-12212,-12126,
-12122,-12118,-12116,-12114,-12112,-11216,-11212,-11116,-11114,-11112,-10325,-10323,-10321,-10315,-10313,-10311,
-10215,-10213,-10211,-5332,-5232,-5222,-5212,-5132,-5122,-5112,-4332,-4232,-4222,-4212,-4132,-4122,-4112,-3334,
-3324,-3322,-3314,-3312,-3303,-3228,-3226,-3224,-3222,-3218,-3216,-3214,-3212,-3203,-3201,-3128,-3126,-3124,
-3122,-3118,-3116,-3114,-3112,-3103,-3101,-2228,-2226,-2224,-2222,-2218,-2216,-2214,-2212,-2203,-2128,-2126,-2124,
-2122,-2118,-2116,-2114,-2112,-2103,-2101,-1218,-1216,-1214,-1212,-1118,-1116,-1114,-1112,-1103,-541,-531,-521,
-511,-431,-421,-411,-327,-325,-323,-321,-317,-315,-313,-311,-217,-215,-213,-211,-16,-15,-14,-13,-12,-11,-6,-5,-4,
-3,-2,-1,0,0,0,0,1,2,3,4,5,6,11,12,13,14,15,16,21,22,111,113,115,117,130,211,213,215,217,221,223,225,227,310,311,
313,315,317,321,323,325,327,331,333,335,337,411,421,431,441,443,511,521,531,541,553,1103,1112,1114,1116,1118,1212,
1214,1216,1218,2101,2103,2112,2114,2116,2118,2122,2124,2126,2128,2203,2212,2214,2216,2218,2222,2224,2226,2228,
3101,3103,3112,3114,3116,3118,3122,3124,3126,3128,3201,3203,3212,3214,3216,3218,3222,3224,3226,3228,3303,3312,
3314,3322,3324,3334,4112,4122,4132,4212,4222,4232,4332,5112,5122,5132,5212,5222,5232,5332,10111,10113,10115,
10211,10213,10215,10221,10223,10225,10311,10313,10315,10321,10323,10325,10331,10333,10335,11112,11114,11116,11212,
11216,12112,12114,12116,12118,12122,12126,12212,12214,12216,12218,12222,12224,12226,13112,13114,13116,13122,13124,
13126,13212,13214,13216,13222,13224,13226,13314,13316,13324,13326,20113,20213,20223,20313,20323,20333,21112,21114,
21212,21214,22112,22114,22122,22124,22212,22214,22222,22224,23112,23114,23122,23124,23126,23212,23214,23222,23224,
23314,23324,30113,30213,30223,30313,30323,31114,31214,32112,32114,32124,32212,32214,32224,33122,33314,33324,41214,
42112,42124,42212,43122,52114,52214,53122,100111,100113,100211,100213,100221,100223,100311,100313,100315,100321,
100323,100325,100331,100333,9000111,9000211,9000221,9010221,9020221,9030221,9030225,9060225,100002110,100002210,
100012110,100012210,1000010020,1000010030,1000020030,1000020040};


static const char* kPnames[FNPART]={
"anti_alpha","anti_He3","anti_triton","anti_deuteron","anti_N(2250)+","anti_N(2250)0","anti_N(2220)+",
"anti_N(2220)0","a0(980)-","k2_star(1980)-","k_star(1410)-","k(1460)-","anti_k2_star(1980)0","anti_k_star(1410)0",
"anti_k(1460)0","rho(1450)-","pi(1300)-","anti_lambda(1810)","anti_N(2090)+","anti_N(2090)0","anti_lambda(1800)",
"anti_N(1710)+","anti_N(1900)+","anti_N(1710)0","anti_N(1900)0","anti_xi(1950)0","anti_xi(1950)-",
"anti_lambda(1670)","anti_delta(1600)++","anti_delta(1600)+","anti_N(1650)+","anti_N(1720)+","anti_delta(1600)0",
"anti_N(1650)0","anti_N(1720)0","anti_delta(1600)-","k_star(1680)-","anti_k_star(1680)0","rho(1700)-",
"anti_xi(1690)0","anti_xi(1690)-","anti_sigma(1940)+","anti_sigma(1750)+","anti_sigma(1940)0","anti_sigma(1750)0",
"anti_lambda(2110)","anti_lambda(1890)","anti_lambda(1600)","anti_sigma(1940)-","anti_sigma(1750)-",
"anti_delta(1920)++","anti_delta(1910)++","anti_delta(1920)+","anti_N(1535)+","anti_N(1700)+","anti_delta(1910)+",
"anti_delta(1920)0","anti_N(1535)0","anti_N(1700)0","anti_delta(1910)0","anti_delta(1920)-","anti_delta(1910)-",
"k1(1400)-","anti_k1(1400)0","a1(1260)-","anti_xi(2030)0","anti_xi(1820)0","anti_xi(2030)-","anti_xi(1820)-",
"anti_sigma(1915)+","anti_sigma(1670)+","anti_sigma(1660)+","anti_sigma(1915)0","anti_sigma(1670)0",
"anti_sigma(1660)0","anti_lambda(1830)","anti_lambda(1690)","anti_lambda(1405)","anti_sigma(1915)-",
"anti_sigma(1670)-","anti_sigma(1660)-","anti_delta(1930)++","anti_delta(1700)++","anti_delta(1900)++",
"anti_N(1990)+","anti_N(1680)+","anti_delta(1700)+","anti_N(1440)+","anti_delta(1930)+","anti_delta(1900)+",
"anti_N(1990)0","anti_N(1680)0","anti_delta(1700)0","anti_N(1440)0","anti_delta(1930)0","anti_delta(1900)0",
"anti_delta(1930)-","anti_delta(1700)-","anti_delta(1900)-","k2(1770)-","k1(1270)-","k0_star(1430)-",
"anti_k2(1770)0","anti_k1(1270)0","anti_k0_star(1430)0","pi2(1670)-","b1(1235)-","a0(1450)-","anti_omega_b-",
"anti_xi_b0","anti_sigma_b+","anti_sigma_b0","anti_xi_b-","anti_lambda_b","anti_sigma_b-","anti_omega_c0",
"anti_xi_c+","anti_sigma_c++","anti_sigma_c+","anti_xi_c0","anti_lambda_c+","anti_sigma_c0","anti_omega-",
"anti_xi(1530)0","anti_xi0","anti_xi(1530)-","anti_xi-","anti_ss1_diquark","anti_sigma(2030)+",
"anti_sigma(1775)+","anti_sigma(1385)+","anti_sigma+","anti_sigma(2030)0","anti_sigma(1775)0","anti_sigma(1385)0",
"anti_sigma0","anti_su1_diquark","anti_su0_diquark","anti_lambda(2100)","anti_lambda(1820)","anti_lambda(1520)",
"anti_lambda","anti_sigma(2030)-","anti_sigma(1775)-","anti_sigma(1385)-","anti_sigma-","anti_sd1_diquark",
"anti_sd0_diquark","anti_delta(1950)++","anti_delta(1905)++","anti_delta++","anti_delta(1620)++",
"anti_delta(1950)+","anti_N(1675)+","anti_delta+","anti_proton","anti_uu1_diquark","anti_N(2190)+",
"anti_delta(1905)+","anti_N(1520)+","anti_delta(1620)+","anti_delta(1950)0","anti_N(1675)0","anti_delta0",
"anti_neutron","anti_ud1_diquark","anti_ud0_diquark","anti_N(2190)0","anti_delta(1905)0","anti_N(1520)0",
"anti_delta(1620)0","anti_delta(1950)-","anti_delta(1905)-","anti_delta-","anti_delta(1620)-","anti_dd1_diquark",
"Bc-","anti_Bs0","B-","anti_B0","Ds-","anti_D0","D-","k3_star(1780)-","k2_star(1430)-","k_star-","kaon-",
"anti_k3_star(1780)0","anti_k2_star(1430)0","anti_k_star0","anti_kaon0","rho3(1690)-","a2(1320)-","rho-","pi-",
"anti_nu_tau","tau+","anti_nu_mu","mu+","anti_nu_e","e+","anti_t_quark","anti_b_quark","anti_c_quark",
"anti_s_quark","anti_u_quark","anti_d_quark","GenericIon","geantino","chargedgeantino","opticalphoton",
"d_quark","u_quark","s_quark","c_quark","b_quark","t_quark","e-","nu_e","mu-","nu_mu","tau-","nu_tau","gluon",
"gamma","pi0","rho0","a2(1320)0","rho3(1690)0","kaon0L","pi+","rho+","a2(1320)+","rho3(1690)+","eta","omega",
"f2(1270)","omega3(1670)","kaon0S","kaon0","k_star0","k2_star(1430)0","k3_star(1780)0","kaon+","k_star+",
"k2_star(1430)+","k3_star(1780)+","eta_prime","phi","f2_prime(1525)","phi3(1850)","D+","D0","Ds+","etac",
"J/psi","B0","B+","Bs0","Bc+","Upsiron","dd1_diquark","delta(1620)-","delta-","delta(1905)-","delta(1950)-",
"delta(1620)0","N(1520)0","delta(1905)0","N(2190)0","ud0_diquark","ud1_diquark","neutron","delta0","N(1675)0",
"delta(1950)0","delta(1620)+","N(1520)+","delta(1905)+","N(2190)+","uu1_diquark","proton","delta+","N(1675)+",
"delta(1950)+","delta(1620)++","delta++","delta(1905)++","delta(1950)++","sd0_diquark","sd1_diquark","sigma-",
"sigma(1385)-","sigma(1775)-","sigma(2030)-","lambda","lambda(1520)","lambda(1820)","lambda(2100)","su0_diquark",
"su1_diquark","sigma0","sigma(1385)0","sigma(1775)0","sigma(2030)0","sigma+","sigma(1385)+","sigma(1775)+",
"sigma(2030)+","ss1_diquark","xi-","xi(1530)-","xi0","xi(1530)0","omega-","sigma_c0","lambda_c+","xi_c0",
"sigma_c+","sigma_c++","xi_c+","omega_c0","sigma_b-","lambda_b","xi_b-","sigma_b0","sigma_b+","xi_b0","omega_b-",
"a0(1450)0","b1(1235)0","pi2(1670)0","a0(1450)+","b1(1235)+","pi2(1670)+","f0(1370)","h1(1170)","eta2(1645)",
"k0_star(1430)0","k1(1270)0","k2(1770)0","k0_star(1430)+","k1(1270)+","k2(1770)+","f0(1710)","h1(1380)",
"eta2(1870)","delta(1900)-","delta(1700)-","delta(1930)-","delta(1900)0","delta(1930)0","N(1440)0","delta(1700)0",
"N(1680)0","N(1990)0","delta(1900)+","delta(1930)+","N(1440)+","delta(1700)+","N(1680)+","N(1990)+",
"delta(1900)++","delta(1700)++","delta(1930)++","sigma(1660)-","sigma(1670)-","sigma(1915)-","lambda(1405)",
"lambda(1690)","lambda(1830)","sigma(1660)0","sigma(1670)0","sigma(1915)0","sigma(1660)+","sigma(1670)+",
"sigma(1915)+","xi(1820)-","xi(2030)-","xi(1820)0","xi(2030)0","a1(1260)0","a1(1260)+","f1(1285)","k1(1400)0",
"k1(1400)+","f1(1420)","delta(1910)-","delta(1920)-","delta(1910)0","N(1700)0","N(1535)0","delta(1920)0",
"delta(1910)+","N(1700)+","N(1535)+","delta(1920)+","delta(1910)++","delta(1920)++","sigma(1750)-","sigma(1940)-",
"lambda(1600)","lambda(1890)","lambda(2110)","sigma(1750)0","sigma(1940)0","sigma(1750)+","sigma(1940)+",
"xi(1690)-","xi(1690)0","rho(1700)0","rho(1700)+","omega(1650)","k_star(1680)0","k_star(1680)+","delta(1600)-",
"N(1720)0","N(1650)0","delta(1600)0","N(1720)+","N(1650)+","delta(1600)+","delta(1600)++","lambda(1670)",
"xi(1950)-","xi(1950)0","N(1900)0","N(1710)0","N(1900)+","N(1710)+","lambda(1800)","N(2090)0","N(2090)+",
"lambda(1810)","pi(1300)0","rho(1450)0","pi(1300)+","rho(1450)+","eta(1295)","omega(1420)","k(1460)0",
"k_star(1410)0","k2_star(1980)0","k(1460)+","k_star(1410)+","k2_star(1980)+","eta(1475)","phi(1680)",
"a0(980)0","a0(980)+","f0(600)","f0(980)","eta(1405)","f0(1500)","f2(1810)","f2(2010)","N(2220)0","N(2220)+",
"N(2250)0","N(2250)+","deuteron","triton","He3","alpha"};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   fNPart(0),
   fNPart30(0),
   fPDG(0),
   fPnames(0),
   fNpReac(FNPREA)
{ 
   memset(fPDGReac,0,fNpReac*sizeof(Short_t));
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
   else return fPrName[ip];
}

//___________________________________________________________________
const Char_t *TPartIndex::ProcNameIndex(Int_t pindex) const {
   if(pindex<0 || pindex>=fNProc) return "Unknown";
   return fPrName[pindex];
}

//___________________________________________________________________
void TPartIndex::SetPartTable(char **names, Int_t *PDG, Int_t np) {
   fNPart = np;
   fNPart30 = 30*fNPart;
   fPnames = new char[fNPart30];
   fPDG = new Int_t[fNPart];
   memset(fPnames,0,fNPart30);
   for(Int_t i=0; i<fNPart; ++i) {
      fPDG[i]=PDG[i];
      strncpy(&fPnames[30*i],names[i],29);
   }
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
	 if(strlen(line)+strlen(&fPnames[30*i])+1>119) {
	    printf("%s\n",line);
	    memset(line,0,120);
	 }
	 strcat(line,&fPnames[30*i]);
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
      R__b.ReadClassBuffer(TPartIndex::Class(),this);
      fgPartIndex = this;
   } else {
      R__b.WriteClassBuffer(TPartIndex::Class(),this);
   }
}




