#ifndef __TNudyEndfRecoPoint__
#define __TNudyEndfRecoPoint__

#include <vector>
#include <fstream>
class TNudyEndfNuPh;
class TNudyEndfFissionYield;
class TNudyEndfEnergy;
class TNudyEndfEnergyAng;
class TNudyEndfAng;
class TNudyEndfFile;
class TNudyEndfList;
class TList;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint> matrixint;
typedef std::vector<rowd> matrixd2;
typedef std::vector<std::vector<rowd>> matrixd3;
typedef std::vector<std::vector<std::vector<rowd>>> matrixd4;
typedef std::vector<std::vector<std::vector<std::vector<rowd>>>> matrixd5;

class TNudyEndfRecoPoint {

public:
  TNudyEndfRecoPoint();
  TNudyEndfRecoPoint(int ielemId, const char *irENDF);
  virtual ~TNudyEndfRecoPoint();
  void GetData(int elemid, const char *irENDF);
  double GetSigmaTotal(int elemid, double energyK);
  double GetSigmaPartial(int elemid, int i, double energyK);
  virtual double GetCos4(int elemid, int mt, double energyK);
  virtual double GetCos64(int elemid, int mt, double energyK);
  virtual int GetCos4Lct(int elemid, int mt);
  virtual double GetEnergy5(int elemid, int mt, double energyK);
  virtual double GetCos6(int elemid, int mt, double energyK);
  virtual double GetEnergy6(int elemid, int mt, double energyK);
  virtual double GetQValue(int elemid, int mt);
  virtual double GetMt4(int elemid, int mt);
  virtual double GetMt5(int elemid, int mt);
  virtual double GetMt6(int elemid, int mt);
  virtual int GetLaw6(int ielemId, int mt);
  virtual int GetZd6(int ielemId, int mt);
  virtual int GetAd6(int ielemId, int mt);
  virtual double GetNuTotal(int elemid, double energyK);
  virtual double GetNuPrompt(int elemid, double energyK);
  virtual double GetNuDelayed(int elemid, double energyK);
  virtual double GetFissHeat(int elemid, double energyK);
  virtual double GetFisYield(int elemid, double energyK);
  virtual double GetLambdaD(int elemid, int time);
  virtual double GetDelayedFraction(int ielemId, int mt, double energyK);
  std::fstream out, outtotal;
  std::string outstring, outstringTotal;
  matrixint MtValues; // MT values for which cross-section/ heating values are given  all elements
protected:
  int elemId;
  const char *rENDF;                      // precision for cross-section reconstruction
  matrixd2 eneUni, sigUniT;               // unionization of energy and total cross-section
  matrixd3 sigUniOfMt;                    // sigma for each reaction afte unionization of energy and all elements
  matrixint energyLocMtId;                // MT wise starting energy for cross-section all elements
  matrixd4 cos4OfMts;                     // cosine and pdf from file 4 for each reaction
  matrixd4 cosPdf4OfMts;                  // cosine and pdf from file 4 for each reaction
  matrixd4 cosCdf4OfMts;                  // cosine and cdf from file 4 for each reaction
  matrixd3 energy4OfMts;                  // incident energy in file 4 for each reaction
  matrixint Mt4Values;                    // MT values for which angular distributions are given in file 4
  matrixint Mt4Lct;                       // CM and Lab flag for angular distributions as given in file 4
  matrixd4 energyOut5OfMts;               // cosine and pdf from file 4 for each reaction
  matrixd4 energyPdf5OfMts;               // cosine and pdf from file 4 for each reaction
  matrixd4 energyCdf5OfMts;               // cosine and cdf from file 4 for each reaction
  matrixd4 cos6OfMts;                     // cosine 6 for each reaction and element
  matrixd4 cosin6Pdf, cosin6Cdf;          // pdf cdf 6 for each reaction and element
  matrixd5 energyOut6OfMts;               // energy from file 6 for each reaction
  matrixd5 energyPdf6OfMts;               // pdf from file 6 for each reaction
  matrixd5 energyCdf6OfMts;               // cdf from file 6 for each reaction
  matrixd3 energy5OfMts;                  // incident energy in file 5 for each reaction
  matrixd3 fraction5OfMts;                // fraction for incident energy in file 5 for each reaction
  matrixint Mt5Values;                    // MT values for which angular distributions are given in file 4
  matrixint Law6;                         // law 6 for angular-energy distributions are given in file 6
  matrixint ZD6, AD6;                     // Z, A of law 6 distributions
  matrixd2 eint, nut;                     // total incident energy and nu,  all elements
  matrixd2 einp, nup;                     // prompt incident energy and nu,  all elements
  matrixd2 eind, nud, lambdaD;            // delayed incident energy and nu,  all elements
  matrixd2 einFissHeat, fissHeat;         // fission incident energy and heat,  all elements
  matrixd2 einfId, qvalue;                // incident energy for fission yield
  matrixd3 zafId, pdfYieldId, cdfYieldId; // za and yield fission
  double AWRI;

private:
  void ReadFile2(TNudyEndfFile *file);
  void ReadFile3(TNudyEndfFile *file);
  void fixupTotal(rowd &x1);

  int flagRead = -1;
  double QValue[999];
  double A, AWR, ABN, QX;         // standard ENDF parameters
  int NR, NP, NE; // standard ENDF parameters for range and interpolation
  matrixint Mt4, Mt5,
      Mt6; // MT values for which angular, energy/ angular-energy distributions are given in file 4, 5, 6
  matrixd2 sigmaOfMts;                      // sigma for each reaction
  matrixd2 sigmaUniOfMts;                   // sigma for each reaction afte unionization of energy
  rowint energyLocationMts;                 // MT wise starting energy for cross-section
  rowint MtNumbers, MtNum4, MtNum5, MtNum6; // MT numbers
  rowd energyMts, sigmaMts, qvaluetemp;                // MT numbers for sigma in file3
  rowd eLinearFile3;
  rowd xLinearFile3;
  rowd eneTemp, sigTemp;       // temporary vectors to store energy and sigma
  TNudyEndfAng *recoAng;
  TNudyEndfEnergy *recoEnergy;
  TNudyEndfEnergyAng *recoEnergyAng;
  TNudyEndfNuPh *recoNuPh;
  TNudyEndfFissionYield *recoFissY;
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif

#ifdef USE_ROOT
  ClassDef(TNudyEndfRecoPoint, 1) // class for an ENDF reconstruction
#endif
};
#endif
