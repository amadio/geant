#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iterator>

#include <Riostream.h>
#include <TTree.h>
#include <TCollection.h>
#include <TIterator.h>
#include <TFile.h>
#include <TMath.h>
#include <TString.h>
#include <TClass.h>
#include <TObject.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TCollection.h>
#include <Rtypes.h>
#include <TCanvas.h>
#include <TNtuple.h>
#include <TTreeReader.h>
#include <TKey.h>
#include <TMacro.h>
#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>

#include "Geant/TNudyDB.h"
#include "Geant/TNudyElementTable.h"
#include "Geant/TNudyENDF.h"
#include "Geant/TNudyEndfCont.h"
#include "Geant/TNudyEndfFile.h"
#include "Geant/TNudyEndfINTG.h"
#include "Geant/TNudyEndfList.h"
#include "Geant/TNudyEndfMat.h"
#include "Geant/TNudyEndfTape.h"
#include "Geant/TNudyEndfRecord.h"
#include "Geant/TNudyEndfSec.h"
#include "Geant/TNudyEndfTab1.h"
#include "Geant/TNudyEndfTab2.h"
#include "Geant/TNudyLibrary.h"
#include "Geant/TNudyManager.h"
#include "Geant/TNudySubLibrary.h"
#include "Geant/TVNudyModel.h"

#define PI 2.0 * asin(1.0)
#define sqrt(x) pow(x,0.5)
#define x2(x) (x * x)
#define x3(x) (x * x2(x))
#define x4(x) (x2(x) * x2(x))
#define x5(x) (x * x4(x))
#define x6(x) (x4(x) * x2(x))
#define x8(x) (x6(x) * x2(x))
#define kconst 2.196771e-3
#define Mn 939.565378e+6    // in eV/c^2   1.00866491588   // in u
#define hcross 6.5821192815e-6   // eV.s  1.054571800e-34 // in J.s
//#define kconst sqrt(2.0 * Mn)/hcross       // used for k and k^prime
#define Fac1(x) (1.0 + x2(x))
#define Fac2(x) (9.0 + 3.0 * x2(x) + x4(x))
#define Fac3(x) 225.0 + 50445.0*x2(x)  - 13500.0*x3(x)  + 1386.0*x4(x)  - 60.0*x5(x)  + x6(x)

class TFile;
class TGeoMatrix;
class TGeoVolume;
class TNudyEndfTape;
class TNudyEndfMat;
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfList;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfINTG;
class TNudyManager;

//****************** make a structure for isotope Parameters
class  resonance {

public:  
  int Z, ZA, ZAI, LFW, NER, LRU, LRF, NRO, NAPS, NLS;
  int LRX;
  double eLo, eHi, SPI, AP, rad_a; // Resonance range Low, High, Target Spin (I), Scattering Radius (AP) Channel radius (a)
  double A, AWR, ABN, AWRI, QX;
  double factor_k;
  double delE;             // this is resolution for Energy
  int JMIN, JMAX;
  int totalPoints ;
  std::vector<int> l;
  std::vector<int> NRS;              // no. of resolved resonances
  std::vector<double> Er;            // resolved resonance energy
  std::vector<double> J;             // associated J
  std::vector<double> GJ;
  std::vector<double> Gamma_r;       // total width = Gamma_n + Gamma_g + Gamma_f
  std::vector<double> Gamma_n;       // neutron scattering width
  std::vector<double> Gamma_g;       // Capture width
  std::vector<double> Gamma_f;       // fission width
  std::vector<double> Gamma_x;       // Gamma_g + Gamma_f
  std::vector<std::vector<double> > ER_mesh;
  // std::vector<double> sigma;
  // std::vector<double> sigmap;
  std::vector<std::vector<int> > JMinMax;
  // std::vector<std::vector<double> > Rho;           // rho(E) and rho'(E)
  //std::vector<std::vector<double> > RhoC;        // rhoCap(E) and rhoCap'(E)

  double K_wnum(double x, int isPrime);
  double Gamma_nrE(double x, int ii, int lval);
  double Gamma_nrEP(double x, int ii, int lval);
  double Gamma_xrE(double x, int ii, int lrx);
  double Gamma_rE(double x, int ii, int lval, int lrx);
  double GetERP(double x, int r, int lVal, int isPrime);
  double GetRho(double x, int isDiff, int lVal);
  double GetRhoC(double x, int isDiff, int lVal);

  double calcPhi(double x, int l);
  double calcPhiP(double x, int l);
  double calcShift(double x, int l);
  double calcShiftP(double x, int l);
  double calcPene(double x, int l);
  double calcPeneP(double x, int l);

  void GetData(const char *rENDF);
  void setParamVector(void);         // X
  void GetSigma(int ii, double x, double &sig, double &sigP);
  void ReadResDat4l(int l1, TNudyEndfList *theList);
  void SLBW();
  void Linearize();  // uniform nergy meshing 
  double GetER(int i) {return this->Er[i]; }
};

void setEnv(void) {
  gSystem->Load("libRIO");
  gSystem->Load("libHist");
  gSystem->Load("libGeom");
  gSystem->Load("libGraf");
  gSystem->Load("libGpad");
  gSystem->Load("libEG");
  gSystem->Load("libMathMore");
  gSystem->Load("libNudy");
  gSystem->Load("/usr/lib64/gcc/x86_64-suse-linux/4.8/libgfortran.so");
}

//************* Convert ENDF to ROOT file function ***********************
void makeRootFile(const char* fENDF, const char* rENDF, Int_t op) {
  TNudyENDF *tn = new TNudyENDF(fENDF, rENDF, "recreate");
  tn->SetLogLev(op);
  tn->Process();
}

double resonance::K_wnum(double x, int isPrime) {
  double k; 
  k = (!isPrime) ? factor_k * sqrt(std::abs(x)) 
    : factor_k * 0.5 * (1.0/sqrt(std::abs(x)));
  return k;
}

double resonance::GetRho(double x, int isDiff, int lVal) {
  if (!NAPS && !isDiff) return K_wnum(x, 0) * rad_a;
  if (!NAPS && isDiff) return K_wnum(x, 1) * rad_a;
  if ( NAPS && !isDiff) return K_wnum(x, 0) * AP;
  if ( NAPS && isDiff) return K_wnum(x, 1) * AP;
}
 
double resonance::GetRhoC(double x, int isDiff, int lVal) {
  if (!isDiff) return K_wnum(x, 0) * AP;
  if (isDiff) return K_wnum(x, 1) * AP;
} 

double resonance::calcPhi(double x, int l) {
  x = GetRhoC(x, 0, l);
  switch (l) {
  case 0:
    return x;
  case 1:
    return (x - atan(x));
  case 2:
    return (x - atan(3.0 * x/(2.0 - x2(x))));
  case 3:
    return (x - atan(x * x2(15.0 - x)/(15.0 - 6.0 * x2(x))));
  }
}

double resonance::calcPhiP(double x, int l) {
  double rhoC = GetRhoC(x, 0, l);
  double rhoCPrime = GetRhoC(x, 1, l);
  x = rhoC;
  switch (l) {
  case 0:
    return rhoCPrime;
  case 1:
    return rhoCPrime - rhoCPrime * (1.0/Fac1(x));
  case 2:
    return rhoCPrime * (x4(x)/Fac2(x));
  case 3:
    return rhoCPrime * ((x6(x) - 60.0 * x5(x) + 1292.0 * x4(x) - 13500.0 * x3(x) + 49500.0 * x2(x) + 900.0 * x - 3150.0) / Fac3(x));
  }
}

double resonance::calcShift(double x, int l) {
  x = GetRho(x, 0, l);
  switch (l) {
  case 0:
    return 0.0;
  case 1:
    return (-1.0 / Fac1(x));
  case 2:
    return (-(18.0 + 3.0 * x2(x)) / Fac2(x));
  case 3:
    return (-(675.0 + 90.0 * x2(x) + 6.0 * x4(x))/(x6(x) + 6.0 * x4(x) + 45.0 * x2(x) + 225.0));
  }
}

double resonance::calcShiftP(double x, int l) {
  double rhoP = GetRho(x, 1, l);
  x = GetRho(x, 0, l);
  switch (l) {
  case 0:
    return 0.0;
  case 1:
    return rhoP * ((2.0 * x)/x2(Fac1(x)));
  case 2:
    return -rhoP * ((6.0 * x * (9.0 + 12.0 * x2(x) + x4(x)))/x2(Fac2(x)));
  case 3:
    return rhoP * ((6.0 * x * (3375.0 + 1800.0 * x2(x) + 765.0 * x4(x) + 60.0 * x6(x) + 2.0 * x6(x) * x2(x)))/x2(225.0 + 45.0 * x2(x) + 6.0 * x4(x) + x6(x)));
  }
}

double resonance::calcPene(double x, int l) {
  x = GetRho(x, 0, l);

  switch (l) {
  case 0:
    return x;
  case 1:
    return (x3(x)/Fac1(x));
  case 2:
    return (x5(x)/Fac2(x));
  case 3:
    return ((x * x6(x))/(225.0 + 45.0 * x2(x) + 6.0 * x4(x) + x6(x)));
  }
}

double resonance::calcPeneP(double x, int l) {
  double rho, rhop;

  rho = GetRho(x, 0, l);
  rhop = GetRho(x, 1, l);
  x = rho;
  switch (l) {
  case 0:
    return rhop;
  case 1:
    return rhop * (((3.0 * x2(x))/Fac1(x) - (2.0 * x4(x))/x2(Fac1(x))));
  case 2:
    return rhop * ((x4(x) * (x4(x) + 9.0 * x2(x) + 45.0))/x2(x4(x) + 3.0 * x2(x) + 9.0));
  case 3:
    return rhop * ((x6(x) * (x6(x) + 18.0 * x4(x) + 225.0 * x2(x) + 1575.0))/x2(x6(x) + 6.0 * x4(x) + 45.0 * x2(x) + 225.0));
  }
}

double resonance::GetERP(double x, int r, int lVal, int isPrime) {
  double result = 0.0;
  double er = Er[r];
 switch (isPrime) {
 case 0:
   result =  (!lVal) ? x : x + (calcShift(std::abs(er), lVal) - calcShift(x, lVal))/(2.0 * calcPene(std::abs(er), lVal)) * Gamma_nrE(std::abs(er), r, lVal);
 case 1:
   result = (!lVal) ? 0.0 : - calcShift(x, lVal)/(2.0 * calcPene(std::abs(er), lVal)) * Gamma_nrE(std::abs(er), r, lVal);
  return result;
 }
}

double resonance::Gamma_nrEP(double x, int ii, int lVal) {
  double er = Er[ii];
  return (calcPeneP(x, lVal) * Gamma_n[ii])/calcPene(std::abs(er), lVal);
}

double resonance::Gamma_nrE(double x, int ii, int lval) {
  double er = Er[ii];
  return  (calcPene(x, lval) * Gamma_n[ii])/calcPene(std::abs(er), lval);
}

double resonance::Gamma_xrE(double x, int ii, int lrx) {
  if (!lrx) {
    return Gamma_r[ii];
  } else {
    return Gamma_r[ii] - (Gamma_n[ii] + Gamma_g[ii] + Gamma_f[ii]);
  }
}

double resonance::Gamma_rE(double x, int ii, int lval, int lrx) {
  return Gamma_nrE(x, ii, lval) + Gamma_xrE(x, ii, lrx) + Gamma_g[ii] + Gamma_f[ii];
}

//******** compute sigma and sigma_prime here ***********
void resonance::GetSigma(int ii, double x, double &sig, double &sigP) {
  std::vector<double> row;
  int nrsl;
  int JMIN, JMAX;
  double temp=0.0;
  double fac1=0.0, fac2 = 0.0;
  double lVal = 0.0;
  double phil, sfil, cfil, sfil2, cfil2, s2fil, c2fil;
  double philp, deno=0.0;
  double sumGJGamma=0.0;
  double ER, ERP, ERPP, GNR, GXR;
  double result = 0.0;
  double resultp = 0.0;
  double fac1p1=0.0, fac1p2=0.0, fac2p1=0.0, fac2p2=0.0;
  double fach, fac2p2a, fac2p2b, fac2p2c, fac2p2d, fac2p2e;
  double sumGJp1 = 0.0, sumGJp2 = 0.0;

  ER = Er[ii];

  // First for Er = a;
  double pibyk2 = PI/x2(K_wnum(x, 0));
  double pibyk3 = PI/x3(K_wnum(x, 0));

  for (int lcnt = 0; lcnt < NLS; lcnt++) {
    lVal = l[lcnt];
    nrsl = NRS[lVal];
    phil = calcPhi(x, lVal);
    philp = calcPhiP(x, lVal);
    sfil = sin(phil);
    cfil = cos(phil);
    sfil2 = x2(sfil);
    cfil2 = x2(cfil);
    s2fil = sin(2.0 * phil);
    c2fil = cos(2.0 * phil);
    fac1 = 4.0 * (2.0 * lVal + 1.0) * sfil2;
    fac1p1 = 4.0 * (2.0 * lVal + 1.0) * 2.0 * sfil * cfil * philp;
    fac1p2 = 4.0 * (2.0 * lVal + 1.0) * (-2.0) * sfil2 * K_wnum(x, 1);
    resultp = pibyk2 * fac1p1 + pibyk3 * fac1p2;

    JMIN = JMinMax[lcnt][0];          JMAX = JMinMax[lcnt][1];
    
    sumGJGamma = 0.0;    sumGJp1 = 0.0;     sumGJp2 = 0.0;
    for (int jVal = JMIN; jVal <= JMAX; jVal++) {
      fac2 = 0.0;  fac2p1 = 0.0;   fac2p2 = 0.0;
      for (int r = 0; r < nrsl; r++) {	
	// Gamma_n,r term
	GNR = Gamma_nrE(x, r, lVal);  // Gamma_nr(E)	
	GXR = Gamma_xrE(x, r, LRX);
	// Er and Er' terms
	ERP = GetERP(x, r,  lVal, 0);
	ERPP = GetERP(x, r, lVal, 1);
	deno = x2(x - ERP) + 0.25 * x2(Gamma_rE(x, r, lVal, LRX));
	fach =  ((x2(GNR) * c2fil) - (2.0 * GNR * GXR * sfil2) + (2.0 * (x - ERP) * GNR * s2fil)/deno);
	fac2 += fach;
	fac2p1 += fach;
	fac2p2a = (philp * GNR * (4.0 * c2fil * (x - ERP) - (2.0 * c2fil * GNR))) / deno;
	fac2p2b = (philp * GNR * (4.0 * cfil * sfil * GXR)) / deno;
	fac2p2c = (2.0 * GNR * s2fil * GNR * (1.0 - ERPP))/ deno;
	fac2p2d = (Gamma_nrEP(x, r, lVal) * (-2.0 * s2fil * (1.0-ERP) + 2.0 * c2fil * GNR - 2.0 * GXR * s2fil)) / deno;
	fac2p2e = ((GNR * (GNR * c2fil) - 2.0 * GXR * s2fil + 2.0 * (x - ERP) * s2fil) * (2.0 * x2(1.0 - ERP) + 0.5 * GNR * Gamma_nrEP(x, r, lVal)))/x2(deno);
	fac2p2 += (fac2p2a - fac2p2b + fac2p2c + fac2p2d - fac2p2e);

      } // END of NRSL r-loop
      sumGJGamma += (GJ[jVal] *  fac2);
      sumGJp1 += GJ[jVal] * fac2p1;
      sumGJp2 += GJ[jVal] * fac2p2;
    }  // END of J jVal-loop
    result += pibyk2 * (fac1 + sumGJGamma); 
    resultp += ((-2.0 * pibyk2 * K_wnum(x, 1)) * sumGJp1 + (pibyk2 * sumGJp2));
  }    // END of NLS l-loop
  sig = result;    sigP = resultp;
  // sigma.push_back(result);
  // sigmap.push_back(resultp);
}

double bisection(double &a, double &b, double &c, double sig1, double sig2, double sig3, double sigP1, double sigP2, double sigP3) {
  double epsa = 1.0e-5;
  double tmp = 0.0;
  double prod = 0.0, prodP = 0.0;
  int i = 1;
  double p0 = 0.0, pi = c;
  if ((std::abs(pi - p0) < epsa) || (std::abs(sig3) < epsa)) return pi;
  while ((sig3 * sig1) > epsa) {
    a = pi; 
    b = b;
    i++;
    p0 = pi;
    pi = 0.5 * (a + b);
    if ((std::abs(pi - p0) < epsa) || (std::abs(sig3) < epsa)) return pi;
  }
  while ((sig3 * sig1) < epsa) {
    a = a;
    b = pi;
    i++;
    pi = 0.5 * (a + b);
    if ((std::abs(pi - p0) < epsa) || (std::abs(sig3) < epsa)) return pi;
  }
  
}
  
//************  Uniform Linearization ***************
void resonance::Linearize() {
  double a = 0.0, b = 0.0, c = 0.0, ER = 0.0;
  double siga = 0.0, sigb = 0.0, sigc = 0.0;
  double sigPa = 0.0, sigPb = 0.0, sigPc = 0.0;
  double bisecE = 0.0;

  bool delESet(false);
  std::vector<double> row;
 
  for (int i = 1; i < (totalPoints - 1); i++) {  // first and last point excluded
    ER = Er[i];     //if (ER < 0.0) continue;
    //if (!delESet) {delE = std::abs(ER * 1.5);   delESet = true;}

    a = ER - delE;    b = ER + delE; c = 0.5 * (a + b);

    GetSigma(i, a, siga, sigPa);     
    GetSigma(i, b, sigb, sigPb);
    GetSigma(i, c, sigc, sigPc);
    if (sigPa >= 0.0 || sigPb < 0.0) bisecE = bisection(a, b, c, siga, sigb, sigc, sigPa, sigPb, sigPc);
    

    std:: cout << i << " ER " << ER << " a " << a << " b " << b << " sgn siga " << sigPa << "  sgn sigb " << sigPb << std::endl;
    

  }
}

//************** Single Level Breit Wigner *********************************
void resonance::SLBW() {
  Linearize();

}


void resonance::ReadResDat4l(int lseries, TNudyEndfList *theList){
  int ii, jj;
  int nrsl = NRS[lseries];
  int nrsl0 = (!lseries) ? 0 : NRS[lseries-1];  // last NRS value
  int lVal = l[lseries];
  int cueMat = 0;

  double gjdeno = 4.0 * SPI + 2.0;
  double r0, r1, rc0, rc1;
  double ErVal;
  double tmp = 0.0, oldE = 0.0;

  delE = 1.0e+5;

  for (ii = 0; ii < nrsl; ii++) { // line 6 onwards data
    cueMat = nrsl0 + ii;
    Er.push_back(theList->GetLIST(ii*6+0)); 
    ErVal = Er[cueMat];
    if (ii) {
      tmp = std::abs(ErVal - oldE);
      delE = (tmp < delE) ? tmp : delE;
    }
    oldE = ErVal;
    J.push_back(theList->GetLIST(ii*6+1));
    Gamma_r.push_back(theList->GetLIST(ii*6+2));
    Gamma_n.push_back(theList->GetLIST(ii*6+3));
    Gamma_g.push_back(theList->GetLIST(ii*6+4));
    Gamma_f.push_back(theList->GetLIST(ii*6+5));	
    GJ.push_back((2.0 * J[cueMat])/gjdeno);
  }
}

//******************** Get Resonant PArameter and cross-section data from ENDF file **********
void resonance::GetData(const char *rENDF) {
  std::vector<int> rowI;

  if (!gGeoManager) gGeoManager = new TGeoManager("rENDF Nudy Manager","");
  TFile *rEND = TFile::Open(rENDF);
  if (!rEND || rEND->IsZombie()) printf("Error: TFile :: Cannot open file %s\n", rENDF);
  TKey *rkey = (TKey*)rEND->GetListOfKeys()->First();
  TNudyEndfTape *rENDFVol = (TNudyEndfTape*)rkey->ReadObj();
  TNudyEndfMat *tMat = 0; 
  TList *mats = (TList*)rENDFVol->GetMats(); 
  
  // I am checking for Char_t* as Root Browser is showing as TString
  //std::cout << " I am checking whether fname and fTitle crashes or not for Char_t *...\n";
  Char_t *theKeyName = (Char_t*)rENDFVol->GetName();
  Char_t *theTitle   = (Char_t*)rENDFVol->GetTitle();
  //std::cout << " Key Name :" << theKeyName << "\nTitle: " << theTitle << std::endl;
  //std::cout << "-------------------------------------------------------\n";
  //std::cout << " I am checking whether fname and fTitle crashes or not for TString ...\n";
  TString strKeyName = rENDFVol->GetName();
  TString strTitle   = rENDFVol->GetTitle();
  //std::cout << " Using TString fName: " << strKeyName << "\nTitle: " << strTitle << std::endl;
  //std::cout << "-------------------------------------------------------\n";
  
  int nmats = mats->GetEntries();
  for (int iMat = 0; iMat < nmats; iMat++) {
    tMat = (TNudyEndfMat*)mats->At(iMat);
    //if (!strcmp(tMat->GetName(),"90-Th-232")) {       //94-Pu-238  90-Th-232
    //if (tMat->GetMAT()==9434) {
    if (tMat->GetMAT()==9434) {
      std::cout << " Found Material " << tMat->GetName() << " < MAT # " << tMat->GetMAT() << ">" << std::endl;
      break;
    }
  }  
  int matNumThis = tMat->GetMAT();  
  TNudyEndfFile *thisFile = (TNudyEndfFile*)tMat->GetFile(2);
  TNudyEndfSec *sec = (TNudyEndfSec*)thisFile->GetSec(151); //GetSection(sectionNumber);
  TList *rec = (TList*)sec->GetRecords();  
  //rec->ls();  // this shows structure inside this TList
  TNudyEndfList *theList = (TNudyEndfList*)rec->First();

  while (theList) {
    if (theList==(rec->At(0))) {  // line 2
      ZAI = theList->GetC1(); ABN = theList->GetC2(); LFW = theList->GetL2(); NER = theList->GetN1();
    }
    if (theList==(rec->At(1))) { // line 3
      eLo = theList->GetC1();  eHi = theList->GetC2();  LRU = theList->GetL1(); 
      LRF = theList->GetL2();  NRO = theList->GetN1();  NAPS = theList->GetN2();
    }
    if (theList==(rec->At(2))) { // line 4
      SPI = theList->GetC1(); AP = theList->GetC2();  NLS = theList->GetN1(); 
      //AP *=1.0e-12;
    }
    NRS.resize(NLS, 0);
    int dataLoc;
    if (theList==(rec->At(3))) { //line 5
      totalPoints = 0;
      for (int lseries = 0; lseries < NLS; lseries++) {  
	//std::cout << "Totalpoints: " << totalPoints << std::endl;  	
	if (lseries) {theList=(TNudyEndfList*)(rec->After(theList));}
	AWRI = theList->GetC1();  
	QX = theList->GetC2();  
	l.push_back(theList->GetL1()); 
	LRX = theList->GetL2();
	totalPoints += theList->GetN2();     // C1 C2 L1 L2 N1 N2
	NRS[lseries] = theList->GetN2();
	A = Mn * AWRI; 
	Z = (int)(ZA - A)/1000.0;  
	rad_a = 0.08 + 0.123 * pow(AWRI, (1./3.));
	factor_k = kconst*(AWRI/(AWRI+1.0));
	
	JMIN = (int)(std::abs(std::abs(SPI - l[lseries])) + 0.5);
	JMAX = (int)(SPI + l[lseries] + 0.5);
	int GJnum = (JMAX - JMIN);	  
	rowI.push_back(JMIN);
	rowI.push_back(JMAX);
	JMinMax.push_back(rowI);
	std::vector<int>().swap(rowI);
	ReadResDat4l(lseries, theList);
      }
      SLBW();
    }
    theList=(TNudyEndfList*)(rec->After(theList));    
  }
  return;
}	     

        
int main() {
  // *************** Define Variables ***********************
  resonance thisMat; 
  bool menu = true;
  int choice;
  Int_t matNum;
  const char* particle = "neutron";
  Int_t zVal, aVal, za;   // zVal = Z, aVal = A, za = 100 * Z + A
  Int_t op = 2;           // Option for Root file creation verbosity
  std::vector<std::vector<double> > CSTable;
    
  // *************** Declare filenames ************************
  const char* fENDF = "e6r1nds8.dat";
  const char* rENDF = "test.root";
  const char* mLib = "mem.dat";
  
  // *************** Initialize some parameters *******************	
  // set some value for test
  //zVal = 90;  aVal = 230; iso = 0; reac = 102; temp=293; awr = 228.0600;
  //zVal = 1;  aVal = 1; iso = 0; reac = 102; temp=293; 

  // ***************   Load Libraries *************************
  setEnv();

  while (menu != false){
    choice = 2;  menu = false;
    if (choice != 2) {
    std::cout << "*******************************\n";
    std::cout << " 1 - Convert ENDF to ROOT file.\n";
    std::cout << " 2 - Process Root file for SLBW. \n";
    
    std::cout << " 9 - Exit.\n";
    std::cout << " Enter your choice and press return: ";

    std::cin >> choice;
    }
    

    switch (choice) {
    case 1:
           // *************** Convert ENDF to ROOT file using NUDY::Federico *********
           makeRootFile(fENDF, rENDF, op);   // <------------WORKING
           break;
    case 2:
           // *************** Get Model for Particular Isotope ********************
           thisMat.GetData(rENDF);
           
	   break;
    case 9:
           std::cout << " Exiting ...........\n";
           menu = false;
           break;
    }
  }

  return 0;
}
