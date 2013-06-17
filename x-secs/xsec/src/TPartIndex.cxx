#include "TPartIndex.h"
#include "TString.h"

ClassImp(TPartIndex)

const char* TPartIndex::fPrName[FNPROC]={"Transport","MultScatt","Ionisation","Decay","inElastic",
			   "Elastic","Capture","Brehms","PairProd","Annihilation",
			   "CoulombScatt","Photoel","Compton","Conversion","Capture",
					"Killer","Total"};
const Short_t TPartIndex::fPCode[FNPROC]={1091,2010,2002,6201,4121,4111,4151,2003,2004,2005,2001,
					  2012,2013,2014,4131,7403,999};

static const Int_t kPDG[FNPART]={
521,-521,511,541,-541,531,411,-411,421,431,-431,0,1000020030,443,12212,12112,2124,1214,22212,22112,32212,32112,2216,
2116,12216,12116,22124,21214,42212,42112,32124,31214,42124,41214,12218,12118,52214,52114,2128,1218,100002210,100002110,
100012210,100012110,553,10211,-10211,10111,9000211,-9000211,9000111,20213,-20213,20113,215,-215,115,1000020040,-511,
-531,-421,-1000020030,-12212,-12112,-2124,-1214,-22212,-22112,-32212,-32112,-2216,-2116,-12216,-12116,-22124,-21214,
-42212,-42112,-32124,-31214,-42124,-41214,-12218,-12118,-52214,-52114,-2128,-1218,-100002210,-100002110,-100012210,
-100012110,-1000020040,-5,-4,-1,-1103,-32214,-32224,-31114,-32114,-2122,-2222,-1112,-1212,-12214,-12224,-11114,-12114,
-12122,-12222,-11112,-11212,-2126,-2226,-1116,-1216,-22122,-22222,-21112,-21212,-22214,-22224,-21114,-22114,-12126,
-12226,-11116,-11216,-2218,-2228,-1118,-2118,-2214,-2224,-1114,-2114,-1000010020,-100311,-10311,-10313,-20313,-10315,
-315,-100315,-317,-100313,-30313,-313,-311,-3122,-13122,-3124,-23122,-33122,-13124,-43122,-53122,-3126,-13126,-23124,
-3128,-23126,-5122,-4122,-2112,-12,-14,-16,-3334,-5332,-4332,-2212,-3,-3101,-3103,-3224,-3114,-3214,-13222,-13112,
-13212,-13224,-13114,-13214,-23222,-23112,-23212,-3226,-3116,-3216,-13226,-13116,-13216,-23224,-23114,-23214,-3228,
-3118,-3218,-3222,-3112,-3212,-5222,-5112,-5212,-4212,-4222,-4112,-3303,-3201,-3203,-6,-1000010030,-2,-2101,-2103,
-2203,-3314,-3324,-23314,-23324,-13314,-13324,-33314,-33324,-13316,-13326,-3312,-3322,-5132,-5232,-4232,-4132,10213,
-10213,10113,5,4,0,1,1103,32214,32224,31114,32114,2122,2222,1112,1212,12214,12224,11114,12114,12122,12222,11112,11212,
2126,2226,1116,1216,22122,22222,21112,21212,22214,22224,21114,22114,12126,12226,11116,11216,2218,2228,1118,2118,2214,
2224,1114,2114,1000010020,-11,11,221,100221,9020221,100331,10225,10335,331,441,10221,9030221,10331,9000221,9010221,
20223,20333,225,9030225,9060225,335,22,0,21,10223,10333,100321,-100321,100311,10321,-10321,10311,10323,-10323,10313,
20323,-20323,20313,10325,-10325,10315,325,-325,315,100325,-100325,100315,327,-327,317,100323,-100323,100313,30323,
-30323,30313,323,-323,313,321,-321,311,130,310,3122,13122,3124,23122,33122,13124,43122,53122,3126,13126,23124,3128,
23126,5122,4122,-13,13,2112,12,14,16,223,100223,30223,3334,227,5332,4332,0,333,100333,337,100211,-100211,100111,211,
-211,111,10215,-10215,10115,2212,100213,-100213,100113,30213,-30213,30113,213,-213,113,217,-217,117,3,3101,3103,3224,
3114,3214,13222,13112,13212,13224,13114,13214,23222,23112,23212,3226,3116,3216,13226,13116,13216,23224,23114,23214,
3228,3118,3218,3222,3112,3212,5222,5112,5212,4212,4222,4112,3303,3201,3203,6,-15,15,1000010030,2,2101,2103,2203,3314,
3324,23314,23324,13314,13324,33314,33324,13316,13326,3312,3322,5132,5232,4232,4132};

static const char* kPnames[FNPART]={
"B+","B-","B0","Bc+","Bc-","Bs0","D+","D-","D0","Ds+","Ds-","GenericIon","He3","J/psi","N(1440)+","N(1440)0",
"N(1520)+","N(1520)0","N(1535)+","N(1535)0","N(1650)+","N(1650)0","N(1675)+","N(1675)0","N(1680)+","N(1680)0",
"N(1700)+","N(1700)0","N(1710)+","N(1710)0","N(1720)+","N(1720)0","N(1900)+","N(1900)0","N(1990)+","N(1990)0",
"N(2090)+","N(2090)0","N(2190)+","N(2190)0","N(2220)+","N(2220)0","N(2250)+","N(2250)0","Upsiron","a0(1450)+",
"a0(1450)-","a0(1450)0","a0(980)+","a0(980)-","a0(980)0","a1(1260)+","a1(1260)-","a1(1260)0","a2(1320)+","a2(1320)-",
"a2(1320)0","alpha","anti_B0","anti_Bs0","anti_D0","anti_He3","anti_N(1440)+","anti_N(1440)0","anti_N(1520)+",
"anti_N(1520)0","anti_N(1535)+","anti_N(1535)0","anti_N(1650)+","anti_N(1650)0","anti_N(1675)+","anti_N(1675)0",
"anti_N(1680)+","anti_N(1680)0","anti_N(1700)+","anti_N(1700)0","anti_N(1710)+","anti_N(1710)0","anti_N(1720)+",
"anti_N(1720)0","anti_N(1900)+","anti_N(1900)0","anti_N(1990)+","anti_N(1990)0","anti_N(2090)+","anti_N(2090)0",
"anti_N(2190)+","anti_N(2190)0","anti_N(2220)+","anti_N(2220)0","anti_N(2250)+","anti_N(2250)0","anti_alpha",
"anti_b_quark","anti_c_quark","anti_d_quark","anti_dd1_diquark","anti_delta(1600)+","anti_delta(1600)++",
"anti_delta(1600)-","anti_delta(1600)0","anti_delta(1620)+","anti_delta(1620)++","anti_delta(1620)-",
"anti_delta(1620)0","anti_delta(1700)+","anti_delta(1700)++","anti_delta(1700)-","anti_delta(1700)0",
"anti_delta(1900)+","anti_delta(1900)++","anti_delta(1900)-","anti_delta(1900)0","anti_delta(1905)+",
"anti_delta(1905)++","anti_delta(1905)-","anti_delta(1905)0","anti_delta(1910)+","anti_delta(1910)++",
"anti_delta(1910)-","anti_delta(1910)0","anti_delta(1920)+","anti_delta(1920)++","anti_delta(1920)-",
"anti_delta(1920)0","anti_delta(1930)+","anti_delta(1930)++","anti_delta(1930)-","anti_delta(1930)0",
"anti_delta(1950)+","anti_delta(1950)++","anti_delta(1950)-","anti_delta(1950)0","anti_delta+","anti_delta++",
"anti_delta-","anti_delta0","anti_deuteron","anti_k(1460)0","anti_k0_star(1430)0","anti_k1(1270)0","anti_k1(1400)0",
"anti_k2(1770)0","anti_k2_star(1430)0","anti_k2_star(1980)0","anti_k3_star(1780)0","anti_k_star(1410)0",
"anti_k_star(1680)0","anti_k_star0","anti_kaon0","anti_lambda","anti_lambda(1405)","anti_lambda(1520)",
"anti_lambda(1600)","anti_lambda(1670)","anti_lambda(1690)","anti_lambda(1800)","anti_lambda(1810)",
"anti_lambda(1820)","anti_lambda(1830)","anti_lambda(1890)","anti_lambda(2100)","anti_lambda(2110)","anti_lambda_b",
"anti_lambda_c+","anti_neutron","anti_nu_e","anti_nu_mu","anti_nu_tau","anti_omega-","anti_omega_b-","anti_omega_c0",
"anti_proton","anti_s_quark","anti_sd0_diquark","anti_sd1_diquark","anti_sigma(1385)+","anti_sigma(1385)-",
"anti_sigma(1385)0","anti_sigma(1660)+","anti_sigma(1660)-","anti_sigma(1660)0","anti_sigma(1670)+",
"anti_sigma(1670)-","anti_sigma(1670)0","anti_sigma(1750)+","anti_sigma(1750)-","anti_sigma(1750)0",
"anti_sigma(1775)+","anti_sigma(1775)-","anti_sigma(1775)0","anti_sigma(1915)+","anti_sigma(1915)-",
"anti_sigma(1915)0","anti_sigma(1940)+","anti_sigma(1940)-","anti_sigma(1940)0","anti_sigma(2030)+",
"anti_sigma(2030)-","anti_sigma(2030)0","anti_sigma+","anti_sigma-","anti_sigma0","anti_sigma_b+","anti_sigma_b-",
"anti_sigma_b0","anti_sigma_c+","anti_sigma_c++","anti_sigma_c0","anti_ss1_diquark","anti_su0_diquark",
"anti_su1_diquark","anti_t_quark","anti_triton","anti_u_quark","anti_ud0_diquark","anti_ud1_diquark",
"anti_uu1_diquark","anti_xi(1530)-","anti_xi(1530)0","anti_xi(1690)-","anti_xi(1690)0","anti_xi(1820)-",
"anti_xi(1820)0","anti_xi(1950)-","anti_xi(1950)0","anti_xi(2030)-","anti_xi(2030)0","anti_xi-","anti_xi0",
"anti_xi_b-","anti_xi_b0","anti_xi_c+","anti_xi_c0","b1(1235)+","b1(1235)-","b1(1235)0","b_quark","c_quark",
"chargedgeantino","d_quark","dd1_diquark","delta(1600)+","delta(1600)++","delta(1600)-","delta(1600)0","delta(1620)+",
"delta(1620)++","delta(1620)-","delta(1620)0","delta(1700)+","delta(1700)++","delta(1700)-","delta(1700)0",
"delta(1900)+","delta(1900)++","delta(1900)-","delta(1900)0","delta(1905)+","delta(1905)++","delta(1905)-",
"delta(1905)0","delta(1910)+","delta(1910)++","delta(1910)-","delta(1910)0","delta(1920)+","delta(1920)++",
"delta(1920)-","delta(1920)0","delta(1930)+","delta(1930)++","delta(1930)-","delta(1930)0","delta(1950)+",
"delta(1950)++","delta(1950)-","delta(1950)0","delta+","delta++","delta-","delta0","deuteron","e+","e-","eta",
"eta(1295)","eta(1405)","eta(1475)","eta2(1645)","eta2(1870)","eta_prime","etac","f0(1370)","f0(1500)","f0(1710)",
"f0(600)","f0(980)","f1(1285)","f1(1420)","f2(1270)","f2(1810)","f2(2010)","f2_prime(1525)","gamma","geantino","gluon",
"h1(1170)","h1(1380)","k(1460)+","k(1460)-","k(1460)0","k0_star(1430)+","k0_star(1430)-","k0_star(1430)0","k1(1270)+",
"k1(1270)-","k1(1270)0","k1(1400)+","k1(1400)-","k1(1400)0","k2(1770)+","k2(1770)-","k2(1770)0","k2_star(1430)+",
"k2_star(1430)-","k2_star(1430)0","k2_star(1980)+","k2_star(1980)-","k2_star(1980)0","k3_star(1780)+","k3_star(1780)-",
"k3_star(1780)0","k_star(1410)+","k_star(1410)-","k_star(1410)0","k_star(1680)+","k_star(1680)-","k_star(1680)0",
"k_star+","k_star-","k_star0","kaon+","kaon-","kaon0","kaon0L","kaon0S","lambda","lambda(1405)","lambda(1520)",
"lambda(1600)","lambda(1670)","lambda(1690)","lambda(1800)","lambda(1810)","lambda(1820)","lambda(1830)",
"lambda(1890)","lambda(2100)","lambda(2110)","lambda_b","lambda_c+","mu+","mu-","neutron","nu_e","nu_mu","nu_tau",
"omega","omega(1420)","omega(1650)","omega-","omega3(1670)","omega_b-","omega_c0","opticalphoton","phi","phi(1680)",
"phi3(1850)","pi(1300)+","pi(1300)-","pi(1300)0","pi+","pi-","pi0","pi2(1670)+","pi2(1670)-","pi2(1670)0","proton",
"rho(1450)+","rho(1450)-","rho(1450)0","rho(1700)+","rho(1700)-","rho(1700)0","rho+","rho-","rho0","rho3(1690)+",
"rho3(1690)-","rho3(1690)0","s_quark","sd0_diquark","sd1_diquark","sigma(1385)+","sigma(1385)-","sigma(1385)0",
"sigma(1660)+","sigma(1660)-","sigma(1660)0","sigma(1670)+","sigma(1670)-","sigma(1670)0","sigma(1750)+",
"sigma(1750)-","sigma(1750)0","sigma(1775)+","sigma(1775)-","sigma(1775)0","sigma(1915)+","sigma(1915)-",
"sigma(1915)0","sigma(1940)+","sigma(1940)-","sigma(1940)0","sigma(2030)+","sigma(2030)-","sigma(2030)0","sigma+",
"sigma-","sigma0","sigma_b+","sigma_b-","sigma_b0","sigma_c+","sigma_c++","sigma_c0","ss1_diquark","su0_diquark",
"su1_diquark","t_quark","tau+","tau-","triton","u_quark","ud0_diquark","ud1_diquark","uu1_diquark","xi(1530)-",
"xi(1530)0","xi(1690)-","xi(1690)0","xi(1820)-","xi(1820)0","xi(1950)-","xi(1950)0","xi(2030)-","xi(2030)0","xi-",
"xi0","xi_b-","xi_b0","xi_c+","xi_c0"};

TPartIndex* TPartIndex::fgPartIndex=0;

//___________________________________________________________________
TPartIndex::TPartIndex():
   fNPart(0),
   fNPart30(0),
   fPDG(0),
   fReacPart(0),
   fPnames(0),
   fNpReac(FNPREA)
{ 
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
   fReacPart = new Int_t[fNPart];
   memset(fPnames,0,fNPart30);
   for(Int_t i=0; i<fNPart; ++i) {
      fPDG[i]=PDG[i];
      strncpy(&fPnames[30*i],names[i],29);
      fReacPart[i]=-1;
   }
}

//______________________________________________________________________________
Bool_t TPartIndex::SetPartReac(const Int_t reacpdg[], Int_t np) {
   if(np != fNpReac) {
      Error("SetPartReac","np is %d and should be %d\n",np,fNpReac);
      exit(1);
      return kFALSE;
   } else {
      for(Int_t i=0; i<np; ++i) 
	 fReacPart[PartIndex(reacpdg[i])]=i;
      return kTRUE;
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




