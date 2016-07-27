#ifndef __TNudyEndfSigma__
#define __TNudyEndfSigma__

#include <vector>
#include <fstream>
class TNudyEndfDoppler;
class TNudyEndfFile;
class TNudyEndfList;
class TList;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

#define PI acos(-1.0)
#define x2(x) (x * x)
#define x3(x) (x * x2(x))
#define x4(x) (x2(x) * x2(x))
#define x5(x) (x * x4(x))
#define x6(x) (x4(x) * x2(x))
#define x8(x) (x6(x) * x2(x))
#define kconst 2.196771e-3
#define Mn 939.565378e+6       // in eV/c^2   1.00866491588   // in u
#define hcross 6.5821192815e-6 // eV.s  1.054571800e-34 // in J.s
#define Fac1(x) (1.0 + x2(x))
#define Fac2(x) (9.0 + 3.0 * x2(x) + x4(x))
#define Fac3(x) 225.0 + 50445.0 * x2(x) - 13500.0 * x3(x) + 1386.0 * x4(x) - 60.0 * x5(x) + x6(x)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint> matrixint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;

class TNudyEndfSigma {

public:
  TNudyEndfSigma();
  TNudyEndfSigma(const char *irENDF, double isigDiff);
  virtual ~TNudyEndfSigma();
  void GetData(const char *irENDF, double isigDiff);
  double SetsigPrecision(double x1) { return sigDiff = x1; }
  std::fstream out, outtotal;
  std::string outstring, outstringTotal;

protected:
private:
  void ReadFile1(TNudyEndfFile *file);
  void ReWriteFile1(TNudyEndfFile *file);
  void ReadFile2(TNudyEndfFile *file);
  void ReadFile3(TNudyEndfFile *file);
  void ReWriteFile3(TNudyEndfFile *file);
  void ReadFile4(TNudyEndfFile *file);
  void ReWriteFile4(TNudyEndfFile *file);
  double recursionLinearFile3(double x1, double x2, double sig1, double sig2, rowd x3, rowd x4);
  void ReadResDat4l(int l1, int mm, TNudyEndfList *theList, TList *rec);
  void recoPlusBroad(int flagNer);
  void GetSigma(int lrfp, double x, double &siga, double &sigb, double &sigc);
  void Linearize(int flagNer);
  void GetSigmaRMP(double x, double &siga, double &sigb, double &sigc);
  void InverseMatrix();
  double backCrsAdler(double x, int l1);
  double calcPhi(double x, int l);
  double calcShift(double x, int l);
  double calcPene(double x, int l);
  double GetERP(double x, int r, int lVal);
  double GetRho(double x, int lVal);
  double GetRhoC(double x, int isDiff, int lVal);
  double Gamma_reduced(double y, int ii, int lval);
  double K_wnum(double x);
  double Gamma_nrE(double x, int ii, int lval);
  double Gamma_xrE(int ii, int lrx);
  double Gamma_rE(double x, int ii, int lval, int lrx);
  void additionalSigma(int LRF, double x);
  double recursionLinear(double x1, double x2, double sig1, double sig2);
  double recursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4, double sig5,
                         double sig6);
  int widthFluctuation(double gnr, double gx, double gg, double gf, int jval);
  double Thinning(rowd &x1, rowd &x2);
  double addFile3Resonance(double &x1, double &x2, rowd &x3, rowd &x4);
  double insertFile3(rowd &x1, rowd &x2);
  double insertFile3High(rowd &x1, rowd &x2);
  void broadSigma(rowd &x1, rowd &x2, rowd &x3);
  void fixupTotal(rowd &x1, rowd &x2);
  double recursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x, std::vector<double> sig);
  double recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  void fillPdf1d();
  void fillPdf2d();
  const char *rENDF;					// precision for cross-section reconstruction
  double sigDiff;					// precision for cross-section reconstruction
  matrixd4 cos4OfMts;        // cosine and pdf from file 4 for each reaction
  matrixd4 cosPdf4OfMts;        // cosine and pdf from file 4 for each reaction
  matrixd4 cosCdf4OfMts;        // cosine and cdf from file 4 for each reaction
  matrixd3 energy4OfMts;       // incident energy in file 4 for each reaction
  matrixint Mt4Values;             // MT values for which angular distributions are given in file 4
  matrixint Mt4Lct;                // CM and Lab flag for angular distributions as given in file 4
  matrixd4 energyOut5OfMts;        // cosine and pdf from file 4 for each reaction
  matrixd4 energyPdf5OfMts;        // cosine and pdf from file 4 for each reaction
  matrixd4 energyCdf5OfMts;        // cosine and cdf from file 4 for each reaction
  matrixd4 cos6OfMts;        // cosine 6 for each reaction and element
  matrixd4 cosin6Pdf, cosin6Cdf; //pdf cdf 6 for each reaction and element
  matrixd5 energyOut6OfMts;        // energy from file 6 for each reaction
  matrixd5 energyPdf6OfMts;        // pdf from file 6 for each reaction
  matrixd5 energyCdf6OfMts;        // cdf from file 6 for each reaction
  matrixd3 energy5OfMts;       // incident energy in file 5 for each reaction
  matrixd3 fraction5OfMts;       // fraction for incident energy in file 5 for each reaction
  matrixint Mt5Values;             // MT values for which angular distributions are given in file 4
  rowd eintFile1,nutFile1,einFile1,nuFile1;
  rowd eindFile1,nudFile1,einphFile1,phFile1;
  rowd einfFile1,heatFile1;
  rowd cnc, nui;
  matrixd2 eint, nut;              // total incident energy and nu,  all elements
  matrixd2 einp, nup;              // prompt incident energy and nu,  all elements
  matrixd2 eind, nud, lambdaD;              // delayed incident energy and nu,  all elements
  matrixd2 einFissHeat, fissHeat;              // fission incident energy and heat,  all elements
  matrixd2 einfId, qvalue;              // incident energy for fission yield
  matrixd3 zafId, pdfYieldId, cdfYieldId;              // za and yield fission 
  double AWRI;
  int Z, ZA, ZAI, LFW, NER, LRU, LRF, NRO, NAPS, NLS, LSSF, NLS2, NJS, INT, NIS,
      intLinLru1 = 0;  // standard ENDF parameters
  int LRX, cueMat = 0; // flag for inelastic reaction, number of J
  int flagRead = -1;
  double QValue[999];
  int dopplerBroad = 0, flagResolve = 0,
      flagUnResolve = 0; // flag for Doppler broadening for thinning, Resolve  and URR parameter exist
  double eLo1 = 0, eLo2 = 0, eHi1 = 0, eHi2 = 0, eLo = 0, eHi = 0; // Resonance energy range Low, High
  double SPI, AP, APL[10], rad_a; // Target Spin (I), Scattering Radius (AP), L-dependent AP, Channel radius (a)
  double A, AWR, ABN, QX;         // standard ENDF parameters
  double factor_k;                // factor for wave vector
  double JMIN, JMAX;              // J values
  double RN, RG, RF, RX;          // standard ENDF parameters
  int totalAdler, crsAdler[4];    // adler adler cross-section
  double R[3][3], S[3][3];        // matrix for RM formalism
  double RI[3][3], SI[3][3];      // matrix for RM formalism
  double MissingJ[5][50], MisGj[5];
  int NJValue[5]; // J values in sorted form
  int NR, NP, NE; // standard ENDF parameters for range and interpolation
  matrixint Mt4, Mt5,
      Mt6; // MT values for which angular, energy/ angular-energy distributions are given in file 4, 5, 6
  rowd energyUni, sigmaUniTotal;            // unionization of energy and total cross-section
  matrixd2 sigmaOfMts;                      // sigma for each reaction
  matrixd2 sigmaUniOfMts;                   // sigma for each reaction afte unionization of energy
  rowint energyLocationMts;                 // MT wise starting energy for cross-section
  rowint MtNumbers, MtNum4, MtNum5, MtNum6; // MT numbers
  rowd sigmaMts, qvaluetemp;                // MT numbers for sigma in file3
  rowd eLinElastic, eLinCapture, eLinFission;
  rowd xLinElastic, xLinCapture, xLinFission;
  rowd eBroadElastic, eBroadCapture, eBroadFission;
  rowd xBroadElastic, xBroadCapture, xBroadFission;
  rowint nbt1, int1;
  rowd eLinearFile3;
  rowd xLinearFile3;
  rowd sigma; 
  rowd ein,cos4,cdf,pdf,lCoef1;
  matrixd2 cos2d,pdf2d,cdf2d,lCoef,ein2d;
  matrixd3 cos3d,pdf3d,cdf3d;
  rowd cosFile4;
  rowd cosPdfFile4;
  rowd cosCdfFile4;
  rowint MtLct;				// LCT numbers
  rowint l;					// l values
  rowint NRS;              			// no. of resolved resonances
  rowint NRJ;              			// no. of URR J
  rowint JSM;              			// URR J
  rowd Er;            			// resolved resonance energy
  rowd J;            			// associated J
  rowd GJ;				// spin multiplication factor
  rowd Gamma_r;       			// total width = Gamma_n + Gamma_g + Gamma_f
  rowd Gamma_n;       			// neutron scattering width
  rowd Gamma_g;       			// Capture width
  rowd Gamma_f;       			// fission width
  rowd Gamma_x;       			// Inelastic width
  rowd Gamma_fa,Gamma_fasq;       	// fission width 1
  rowd Gamma_fb,Gamma_fbsq;       	// fission width 2
  rowd at1;       			// 1 background constant (Reich-Moore)
  rowd at2;       			// 2 background constant (Reich-Moore)
  rowd at3;       			// 3 background constant (Reich-Moore)
  rowd at4;       			// 4 background constant (Reich-Moore)
  rowd bt1;       			// 5 background constant (Reich-Moore)
  rowd bt2;       			// 6 background constant (Reich-Moore)
  rowd det1;       			// 1 resonance energy (Reich-Moore)
  rowd dwt1;       			// 2 half width (Reich-Moore)
  rowd grt1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd git1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd def1;       			// 5 background constant (Reich-Moore)
  rowd dwf1;       			// 6 background constant (Reich-Moore)
  rowd grf1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd gif1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd dec1;       			// 5 background constant (Reich-Moore)
  rowd dwc1;       			// 6 background constant (Reich-Moore)
  rowd grc1;       			// 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd gic1;       			// 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd amux, amun, amug, amuf;		// standard ENDF parameters
  rowd Es;				// energy URR
  rowd D, GX, GNO, GG, GF;		// URR parameters
  rowd PhiEr,ShiftEr;			// penetration and shift factors
  rowd eneTemp,sigTemp;			// temporary vectors to store energy and sigma
  TNudyEndfDoppler *doppler;
  TNudyEndfAng *recoAng;
  TNudyEndfEnergy *recoEnergy;
  TNudyEndfEnergyAng *recoEnergyAng;
  TNudyEndfNuPh *recoNuPh;
  TNudyEndfFissionYield *recoFissY;
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyEndfSigma, 1) // class for an ENDF reconstruction
#endif
};
#endif
