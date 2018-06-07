//===-- Nudy/TNudyEndfSigma.h - Instruction class definition -------*- C++ -*-===//
//
//                     The Project Nudy
//===----------------------------------------------------------------------===//
///
/// \file This class processes the neutron cross-section from endf file 2 and 3
///  makes it linearly interpolable, constructs resonance cross-section
/// \brief The class has declaration of  TNudyEndfDoppler which does doppler broadening
/// TNudyEndfAng processes energy integrated angular distribution for neutron
/// TNudyEndfEnergy processes angle integrated energy distribution for neutron
/// TNudyEndfEnergyAng processes energy and angle corelated distribution
/// TNudyEndfNuPh processes fission neutron multiplicity, photon energy distribution
/// from file 1 of the endf format
/// TNudyEndfFissionYield processes fission yield for the fissile material
/// TNudyEndfPhYield processes secondary photon production multiplicity
/// TNudyEndfPhProd production cross-section for the secondary photons
/// TNudyEndfPhAng processes energy integrated angular distribution for photon
/// TNudyEndfPhEnergy processes angle integrated energy distribution for photon
/// TNudyENDF used to store the data in the list format
/// TNudyEndfFile reads the filewise data from the endf format
/// TNudyEndfList reads the list indexed data from the endf format
/// \class TNudyEndfSigma
/// \author H. Kumawat
/// \date July 2016
//===----------------------------------------------------------------------===//
#ifndef __TNudyEndfSigma__
#define __TNudyEndfSigma__

#include <vector>
#include <fstream>

namespace NudyPhysics {
class TNudyEndfDoppler;
class TNudyEndfAng;
class TNudyEndfEnergy;
class TNudyEndfEnergyAng;
class TNudyEndfPhYield;
class TNudyEndfPhProd;
class TNudyEndfPhAng;
class TNudyEndfPhEnergy;
}

namespace Nudy {
class TNudyENDF;
class TNudyEndfFile;
class TNudyEndfList;
class TNudyEndfTab1;
class TNudyEndfTab2;
}

class TList;
class TIter;

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

#define PI acos(-1.0)
#define x2(x) (x * x)
#define x3(x) (x * x2(x))
#define x4(x) (x2(x) * x2(x))
#define x5(x) (x * x4(x))
#define x6(x) (x4(x) * x2(x))
#define x8(x) (x6(x) * x2(x))
#define kconst 2.196771e-3
#define Mn 939.565378e+6
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

namespace NudyPhysics {

class TNudyEndfSigma {

public:
  TNudyEndfSigma();
  /**
   * @param[in] irENDF      ENDF Root file name, this file is generated from the TNudyENDF class.
   * @param[in] isigDiff    Tolerance for the linearization modules in reconstruction/doppler etc.
   */
  TNudyEndfSigma(const char *irENDF, double isigDiff);
  /// \brief constructor to pass root input data file and tolerance  to reconstruct the resonance cross-section
  virtual ~TNudyEndfSigma();
  /// \brief destructor
  void GetData(const char *irENDF, double isigDiff);
  /// \brief this the main function to interact from outside of this class to get the cross-sections
  void SetsigPrecision(double x1) { fSigDiff = x1; }
  /// \brief setting error tolerance for linearly interpolable cross-sections
  void SetInitTempDop(double t1) { fDoppTemp1 = t1; }
  /// \brief set initial temperature for the doppler broadening
  void SetOutTempDop(double t2) { fDoppTemp2 = t2; }
  /// \brief set final temperature for the doppler broadening
  double GetsigPrecision() const { return fSigDiff; }
  /// \brief get precision for the linearly interpolable cross-sections which was set by SetsigPrecision function
  double GetOutTempDop() const { return fDoppTemp2; }
  /// \brief get final temperature for the doppler broadening which was set by SetOutTempDop function
  void SetPreProcess(int x1) { fPrepro = x1; }
  /// \brief set x1 = 0 for basic endf file and x1 = 1 if file is processed by Prepro and written in ENDF format
  virtual void ProcessTab1(Nudy::TNudyEndfTab1 *tab1,int &NR,int &NP,rowint &fnbt, rowint &fint,rowd &x1, rowd &x2);
  /// \brief process tab1 entry
  virtual void ProcessTab2(Nudy::TNudyEndfTab2 *tab2, int &NR,int &NP,rowint &fnbt, rowint &fint); 
  /// \brief process tab2 entry 
  virtual void ModifyTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, double &x);
  /// \brief modify Tab1 record after linearization
  virtual void CreateTab2(Nudy::TNudyEndfTab2 *secTab2, int &NE);
  /// \brief create tab2 record in file5 (MF = 5) to make LF=1 structure
  virtual void CreateTab1(Nudy::TNudyEndfTab1 *secTab1, rowd &x1, rowd &x2, int mtf[], int law, double &x);
  /// \brief create tab1 record in file5 (MF = 5) to make LF=1 structure
  std::fstream out, outtotal;
  std::string outstring, outstringTotal;
  bool GetIsFiss() { return IsFission; }
  
protected:
  matrixd4 fCos4OfMts;
  /// \brief 4-D vector: cosine in file 4 for each given reaction and element
  matrixd4 fCosPdf4OfMts;
  /// \brief cosine pdf in file 4 for each given reaction and element
  matrixd4 fCosCdf4OfMts;
  /// \brief cosine cdf from file 4 for each given reaction and element
  matrixd3 fEnergy4OfMts;
  /// \brief incident energy in file 4 for which cosine angles are given
  matrixint fMt4Values;
  /// \brief MT values for which angular distributions are given in file 4
  matrixint fMt4Lct;
  /// CM and Lab flag for angular distributions as given in file 4
  matrixd4 fEnergyOut5OfMts;
  /// \brief secondary energy for neutron from file 5 for each given reaction and element
  matrixd4 fEnergyPdf5OfMts;
  /// \brief secondary energy pdf for neutron from file 5 for each given reaction and element
  matrixd4 fEnergyCdf5OfMts;
  /// \brief secondary energy cdf for neutron from file 5 for each given reaction and element
  matrixd3 fEnergy5OfMts;
  /// \brief incident neutron energy in file 5 for given reaction and element for which energy distribution is given
  matrixd3 fFraction5OfMts;
  /// \brief delayed neutron fraction for incident energy for given reaction and element
  matrixint fMt5Values;
  /// \brief MT values for which energy distributions are given in file 5
  double AWRI;
  /// \brief mass per neutron, standard ENDF parameters
  matrixd4 fCos6OfMts;
  /// \brief incident cosine from file 6 for given reaction and element
  matrixd4 fCosin6Pdf, fCosin6Cdf;
  /// \brief given cosine pdf  and cdf in file 6 for given reaction and element
  matrixd5 fEnergyOut6OfMts;
  /// \brief energy from file 6 for each reaction and element
  matrixd5 fEnergyPdf6OfMts;
  /// \brief outgoing neutron's pdf from file 6 for given reaction and element
  matrixd5 fEnergyCdf6OfMts;
  /// \brief outgoing neutron's cdf from file 6 for given reaction and element
  matrixint fLaw6;
  /// \brief law 6 for angular-energy distributions are given in file 6
  matrixint fZD6, fAD6;
  /// \brief Z, A of law 6 distributions
  int mtf[3];
  /// \brief MAT, MT, MF parameters
  
private:
  void ReadFile1(Nudy::TNudyEndfFile *file);
  /// \brief reading of the file1 (MF = 1) from root data file that is written by TNudyENDF class
  void ReWriteFile1(Nudy::TNudyEndfFile *file);
  /// \brief rewrite data after interpolating in the linear format for file1 (MF = 1) in root data file
  void ReadFile2(Nudy::TNudyEndfFile *file);
  /// \brief reading of the file2 (MF = 2) from root data file (resonance parameters)
  void ReadFile3(Nudy::TNudyEndfFile *file);
  /// \brief reading of the file3 (MF = 3) from root data file (cross-sections other than resonance cross-sections)
  void ReWriteFile3(Nudy::TNudyEndfFile *file);
  /// \brief rewrite all linearized cross-section data of file 2+3 (MF = 2, 3) into root data file 3
  double RecursionLinearFile3(double x1, double x2, double sig1, double sig2, rowd x3, rowd x4);
  /// \brief linearization program for the cross-section data
  void RecoPlusBroad(int flagNer);
  /// \brief call for reconstruction extra points near thermal region and urr region
  void GetSigma(int lrfp, double x, double &siga, double &sigb, double &sigc);
  /// \brief getting reconstructed cross-section based on different type of formulation
  void Linearize(int flagNer);
  /// \brief call to reconstruct cross-sections for different angular momentum values
  void GetSigmaRMP(double x, double &siga, double &sigb, double &sigc);
  /// \brief getting reconstructed cross-section based on RM formulation
  void InverseMatrix();
  /// \brief matrix inversion
  double BackCrsAdler(double x, int l1);
  /// \brief background cross-section for adler adler formulation
  double CalcPhi(double x, int l);
  /// \brief angle phi claculation for the SLBW and MLBW cross-section calculation
  double CalcShift(double x, int l);
  /// \brief shift in angle phi
  double CalcPene(double x, int l);
  /// \brief penetration factor for different l values
  double GetERP(double x, int r, int lVal);
  /// \brief Shifted resonance energy after phi and penetration factor
  double GetRho(double x, int lVal);
  /// \brief neutron wave number times channel radius parameter \cite ENDF manual
  double GetRhoC(double x, int isDiff, int lVal);
  /// \brief neutron wave number times energy dependent radius parameter
  double Gamma_reduced(double y, int ii, int lval);
  /// \brief reduced gamma widths
  double K_wnum(double x);
  /// \brief neutron wave number
  double Gamma_nrE(double x, int ii, int lval);
  /// \brief neutron width calculation
  double Gamma_xrE(int ii, int lrx);
  /// \brief extra channel width calculation
  double Gamma_rE(double x, int ii, int lval, int lrx);
  /// \brief Total width calculation
  void AdditionalSigma(int LRF, double x);
  /// \brief adding additional cross-section points for linearization to work
  double RecursionLinear(double x1, double x2, double sig1, double sig2);
  /// \brief Recursive linear cross-section between two energy points
  double RecursionLinear(double x1, double x2, double sig1, double sig2, double sig3, double sig4, double sig5,
                         double sig6);
  /// \brief Recursive linear cross-section between two energy points
  int WidthFluctuation(double gnr, double gx, double gg, double gf, int jval);
  /// \brief width fluction parameter for the urr region
  double Thinning(rowd &x1, rowd &x2);
  /// \brief eleminating cross-section points if it is more than linear interpolable and within tolerance
  double AddFile3Resonance(double &x1, double &x2, rowd &x3, rowd &x4);
  /// \brief adding file 3 cross-section data in URR region LSSF=0
  double InsertFile3(rowd &x1, rowd &x2);
  /// \brief adding file 3 cross-section data in URR region LSSF=1
  double InsertFile3High(rowd &x1, rowd &x2);
  /// \brief adding file 3 cross-section data after RR+URR region (high energy region)
  int BinarySearch(double x1, rowd &x2);
  /// \brief binary search algo
  void AddSecFile3(Nudy::TNudyEndfFile *file, double a, double b, int MF2Add, rowd &x1, rowd &x2);
  /// \brief adding integral charge particle cross-sections from MT =103-107 if it does not exist
  /// but exist in terms of individual excited state cross-section
  void BroadSigma(rowd &x1, rowd &x2, rowd &x3);
  /// \brief call to doppler broadening
  void FixupTotal(Nudy::TNudyEndfFile *file, rowd &x1, rowd &x2);
  /// \brief getting total cross-section after doppler broadening
  void DopplerAll();
  /// \brief call to doppler broadening all corss-sections
  double RecursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x,
                             std::vector<double> sig);
  /// \brief linearization of photon cross-sections given in the file1
  double RecursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the legendre polynomial
  double RecursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the probability distribution
  void FillPdf1D();
  /// \brief filling 1 dimentional pdf for anglular distribution
  void FillPdf2D();
  /// \brief filling 2 dimentional pdf for anglular distribution
  void SetIsFission(bool FissKey) { IsFission = FissKey; }
  /// \brief setting if fission takes place
  double RecursionLinearLeg1D(double x1, double x2, double pdf1, double pdf2);
  /// \brief linearization of angular distribution if it is given in terms if the legendre polynomial
  void ConvertFile4ToFile6(Nudy::TNudyEndfFile *file);
  /// \brief function to convert file 4 (angular distribution) to file 6
  void ConvertFile5ToFile6(Nudy::TNudyEndfFile *file);
  /// \brief function to convert file 5 (energy distribution) to file 6
  void ConvertFile14ToFile6(Nudy::TNudyEndfFile *file);
  /// \brief function to convert file 14 (photon angular distribution) to file 6
  void ConvertFile15ToFile6(Nudy::TNudyEndfFile *file);
  /// \brief function to convert file 15 (photon energy distribution) to file 6
  
  const char *rENDF;
  /// \brief Name of the endf to root data file
  double fDoppTemp1, fDoppTemp2;
  /// \brief initial and final temperature for the doppler broadening
  double fSigDiff;
  /// \brief precision/tolerance for cross-section reconstruction while linearization from true values
  matrixd4 fCos4OfMtsPhoton;
  /// \brief 4-D vector: photon cosine in file 4 for each given reaction and element
  matrixd4 fCosPdf4OfMtsPhoton;
  /// \brief cosine photon pdf in file 4 for each given reaction and element
  matrixd4 fCosCdf4OfMtsPhoton;
  /// \brief cosine photon cdf from file 4 for each given reaction and element
  matrixd3 fEnergy4OfMtsPhoton;
  /// \brief incident neutron energy in file 4 for which cosine angles are given
  matrixint fMt4ValuesPhoton;
  /// \brief MT values for which angular distributions of photons are given in file 4
  matrixint fMt4LctPhoton;
  /// CM and Lab flag for angular distributions of photons as given in file 4
  rowd fEintFile1, fNutFile1, fEinFile1, fNuFile1;
  /// \brief energy, total fission neutron multiplicity, energy, prompt fission neutron multiplicity
  rowd fEindFile1, fNudFile1, fEinPhFile1, fPhFile1;
  /// \brief energy, delayed fission neutron multiplicity, energy, prompt fission neutron multiplicity
  rowd fEinfFile1, fHeatFile1;
  /// \brief energy, fission heat
  rowd fCnc, fNui;
  /// \brief coefficients for getting neutron multiplicity \cite ENDF manual
  matrixd2 fEint, fNut;
  /// \brief incident energy and total nu,  all elements
  matrixd2 fEinp, fNup;
  /// \brief prompt incident energy and nu,  all elements
  matrixd2 fEind, fNud, fLambdaD;
  /// \brief delayed incident energy and nu and delay lambda for all elements
  matrixd2 fEinFissHeat, fFissHeat;
  /// \brief fission incident energy and heat,  all elements
  matrixd2 fEinfId, fQvalue;
  /// \brief incident energy for fission yield and qvalues for the reactions
  matrixd3 fZafId, fPdfYieldId, fCdfYieldId;
  /// \brief za, pdf and cdf for fssion yield
  int Z, ZA, ZAI, LFW, NER, LRU, LRF, NRO, NAPS, NLS, fLSSF, NLS2, NJS, INT, NIS, intLinLru1 = 0;
  /// \brief standard ENDF parameters
  bool IsFission;
  /// \brief Abhijit:: to simply know whether fit for fission
  int fLrx, fCueMat = 0;
  /// \brief flag for inelastic reaction, number of J
  int fFlagRead = -1;
  /// \brief flag for reading charge particle production cross-sections
  double fQValue[999];
  /// \brief q-value for all the reactions
  int fDopplerBroad = 0, fFlagResolve = 0, fFlagUnResolve = 0;
  /// \brief flag for thinning before doppler broadening, flas if, Resolve and URR parameter exist
  double fElo1 = 0, fElo2 = 0, fEhi1 = 0, fEhi2 = 0, fElo = 0, fEhi = 0;
  /// \brief Resonance energy range Low, High energy boundaries for RR and URR regions
  double fSpi, fAp, fApl[10], fRad_a;
  /// \brief Target Spin (I), Scattering Radius (AP), L-dependent AP, Channel radius (a)
  double fA, fAWR, fABN, fQX;
  /// \brief mass number mass in units of neutron mass, abundance, q-value from ENDF parameters
  double fFactor_k;
  /// \brief factor for wave vector
  double fJMin, fJMax;
  /// \brief min and max J values
  double fRN, fRG, fRF, fRX;
  /// \brief standard ENDF parameters for URR for neutron, gamma, fission and extra channel
  int totalAdler, crsAdler[4];
  /// \brief adler adler cross-section
  double fR[3][3], fS[3][3];
  /// \brief matrix for RM formalism
  double fRI[3][3], fSI[3][3];
  /// \brief matrix for RM formalism
  double fMissingJ[5][50], fMisGj[5];
  /// \brief J vlaues for calculation of the missing J value factor, spin factor
  int fNJValue[5];
  /// \brief J values in sorted form
  int fNR, fNP, fNE, fMAT, fMT, fMF;
  /// \brief standard ENDF parameters for range and interpolation
  matrixint fMt4, fMt5, fMt6;
  /// \brief MT values for which angular, energy/ angular-energy distributions are given in file 4, 5, 6
  rowd energyUni, sigmaUniTotal;
  /// \brief unionization of energy and total cross-section
  matrixd2 fSigmaOfMts;
  /// \brief sigma for each reaction
  matrixd2 fSigmaOfMtsDop;
  /// \brief sigma for each reaction after doppler
  matrixd2 fSigmaUniOfMts;
  /// \brief sigma for each reaction afte unionization of energy
  rowint fEnergyLocationMts;
  /// \brief MT wise starting energy for cross-section
  rowint MtNumbers, MtNum4, MtNum5, MtNum6;
  /// \brief MT numbers for different distributions in file 4, 5, 6
  rowint MtNumSig4Photon, MtNumAng4Photon, MtNumEng4Photon;
  /// \brief MT numbers for photon cross-section, angle, energy
  rowd fSigmaMts, fQvalueTemp;
  /// \brief MT numbers for sigma in file3 and temp. varaible for q-value
  rowd fELinElastic, fELinCapture, fELinFission;
  /// \brief linear energy points for elastic, capture and fission cross-section
  rowd fXLinElastic, fXLinCapture, fXLinFission;
  /// \brief linear cross-section points for elastic, capture and fission
  rowd fXBroadElastic, fXBroadCapture, fXBroadFission;
  /// \brief linear cross-section points after doppler broadening for elastic, capture and fission
  rowd fE1, fP1, fE2, fP2;
  /// \brief standard file 5 parameters from ENDF manual
  rowint fNbt1, fInt1;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  rowint fNbt2, fInt2;
  /// \brief standard ENDF interpolation parameter \cite ENDF Manual
  int fNr2, fNp2;
  /// \brief standard ENDF parameters for no .of regions and points for interpolation
  rowd fELinearFile3, fXLinearFile3;
  /// \brief linear cross-section points for file 3 data
  rowd fSigma;
  /// \brief temp. cross-section variable for doppler broadening
  rowd fEin, fCos4, fCdf, fPdf, fLegendCoef1;
  /// \brief Temp. variables for cosine, cdf, cdf and legendre coefficient for
  /// angular distribution
  matrixd2 fCos2D, fPdf2D, fCdf2D, fLegendCoef, fEin2D;
  /// \brief Temp. variables for cosine, pdf, cdf and legendre coefficient and
  /// energy for angular distribution
  matrixd3 fCos3D, fPdf3D, fCdf3D;
  /// \brief Temp. variables for cosine, pdf, cdf for angular distribution
  rowd fCosFile4, fCosPdfFile4, fCosCdfFile4;
  /// \brief cosine, pdf and cdf for angular distribution
  rowint fMtLct;
  /// \brief LCT number flag for lab or center of momentum system
  rowint fMtNumbers;
  /// \brief temp. MT numbers
  rowint l; // l values
  rowint fNRS;
  /// \brief no. of resolved resonances
  rowint fNRJ;
  /// \brief no. of URR J values
  rowint fJSM;
  /// \brief URR J values
  rowd fEr;
  /// \brief resolved resonance energy
  rowd J;
  /// \brief associated J
  rowd fGJ;
  /// \brief spin multiplication factor
  rowd fGamma_r;
  /// \brief total width = Gamma_n + Gamma_g + Gamma_f
  rowd fGamma_n;
  /// \brief neutron scattering width
  rowd fGamma_g;
  /// \brief Capture width
  rowd fGamma_f;
  /// \brief fission width
  rowd fGamma_x;
  /// \brief Inelastic width
  rowd fGamma_fa, fGamma_fasq;
  /// \brief fission width 1 for RM formulation \cite ENDF Manual
  rowd fGamma_fb, fGamma_fbsq;
  /// \brief fission width 2 for RM formulation \cite ENDF Manual
  rowd fAt1;
  /// \brief 1st background constant (Reich-Moore)
  rowd fAt2;
  /// \brief 2 background constant (Reich-Moore)
  rowd fAt3;
  // \brief 3 background constant (Reich-Moore)
  rowd fAt4;
  /// \brief 4 background constant (Reich-Moore)
  rowd fBt1;
  /// \brief 5 background constant (Reich-Moore)
  rowd fBt2;
  /// \brief 6 background constant (Reich-Moore)
  rowd fDet1;
  /// \brief 1 resonance energy (Reich-Moore)
  rowd fDwt1;
  /// \brief 2 half width (Reich-Moore)
  rowd fGrt1;
  /// \brief 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd fGit1;
  /// \brief 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd fDef1;
  /// \brief 5 background constant (Reich-Moore)
  rowd fDwf1;
  /// \brief 6 background constant (Reich-Moore)
  rowd fGrf1;
  /// \brief 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd fGif1;
  /// \brief 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd fDec1;
  /// \brief 5 background constant (Reich-Moore)
  rowd fDwc1;
  /// \brief 6 background constant (Reich-Moore)
  rowd fGrc1;
  /// \brief 3 symmetrical cross-section parameter G (Reich-Moore)
  rowd fGic1;
  /// \brief 4 Asymmetrical total cross section parameter, HTr (Reich-Moore)
  rowd fAmux, fAmun, fAmug, fAmuf;
  /// \brief standard ENDF parameters for degrees of freedom in URR for inelastic, elastic, cap., and fiss.
  rowd fEs;
  /// \brief centroid energy URR
  rowd fD, fGX, fGNO, fGG, fGF;
  /// \brief URR parameters \cite ENDF Manual
  rowd fPhiEr, fShiftEr;
  /// \brief penetration and shift factors
  rowd fEneTemp, fSigTemp;
  /// \brief temporary vectors to store energy and sigma
  rowd fEneUniP, fSigUniP;
  /// \brief unionization of energy and total cross-section for n,p
  rowd fEneUniD, fSigUniD;
  /// \brief unionization of energy and total cross-section for n,d
  rowd fEneUniT, fSigUniT;
  /// \brief unionization of energy and total cross-section for n,t
  rowd fEneUniHe3, fSigUniHe3;
  /// \brief unionization of energy and total cross-section for n,He3
  rowd fEneUniHe4, fSigUniHe4;
  /// \brief unionization of energy and total cross-section for n,He4
  rowint fEneLocP, fEneLocD, fEneLocT, fEneLocHe3, fEneLocHe4;
  /// \brief location of the first energy grid point
  matrixd2 fSigUniOfP, fSigUniOfD, fSigUniOfT, fSigUniOfHe3, fSigUniOfHe4;
  /// \brief cross-section at union energy grids for charge particles
  rowd fEneUniPAll, fSigUniPAll;
  /// \brief unionization of energy and total cross-section for n,p + n,pX
  rowd fEneUniDAll, fSigUniDAll;
  /// \brief unionization of energy and total cross-section for n,d + n,dX
  rowd fEneUniTAll, fSigUniTAll;
  /// \brief unionization of energy and total cross-section for n,t + n,tX
  rowd fEneUniHe3All, fSigUniHe3All;
  /// \brief unionization of energy and total cross-section for n,He3 + n,He3X
  rowd fEneUniHe4All, fSigUniHe4All;
  /// \brief unionization of energy and total cross-section for n,He4 + n,He4X
  rowint fEneLocPAll, fEneLocDAll, fEneLocTAll, fEneLocHe3All, fEneLocHe4All;
  /// \brief location of the first energy grid point
  matrixd2 fSigUniOfPAll, fSigUniOfDAll, fSigUniOfTAll, fSigUniOfHe3All, fSigUniOfHe4All;
  /// \brief cross-section at union energy grids for charge particles for all the excited states
  /// and other than charge particle is also emitted
  int fMTChargeFlag[10];
  /// \brief flag if charge particle production is added in MT = 103-107
  int fPrepro = 0;
  /// \brief flag =0 if endf file is processed =1 if prepro processed file is processed (for comparision)
  rowint fNSecNeutron, fNSecPhoton;
  /// \brief counter for neutron, photon producing reactions
  rowint fNReacNeutron, fNDelayFamily;
  /// \brief number of neutron reactions, delayed neutron families
  int fMloop;
  /// \brief termination of recursive loop
  rowd fE4, fP4;
  /// \brief energy, probability temporary variables to fill tab record
  NudyPhysics::TNudyEndfDoppler *fDoppler;
  NudyPhysics::TNudyEndfAng *fRecoAng;
  NudyPhysics::TNudyEndfEnergy *fRecoEnergy;
  NudyPhysics::TNudyEndfEnergyAng *fRecoEnergyAng;
  NudyPhysics::TNudyEndfPhYield *fRecoPhYield;
  NudyPhysics::TNudyEndfPhProd *fRecoPhProd;
  NudyPhysics::TNudyEndfPhAng *fRecoPhAng;
  NudyPhysics::TNudyEndfPhEnergy *fRecoPhEnergy;
  Nudy::TNudyENDF *fPendf;

#ifdef USE_ROOT
  ClassDef(TNudyEndfSigma, 1) // class for an ENDF reconstruction
#endif
};

} // namespace
#endif
