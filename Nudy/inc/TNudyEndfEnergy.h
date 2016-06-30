#ifndef TNudyEndfEnergy_H
#define TNudyEndfEnergy_H

#include "TNudyEndfRecoPoint.h"

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd > matrixd2;
typedef std::vector<std::vector<rowd > > matrixd3;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom;
#endif

class  TNudyEndfEnergy : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergy ();
  TNudyEndfEnergy (TNudyEndfFile *file);
  virtual double GetEnergy5(int elemid, int mt, double energyK);
  virtual ~TNudyEndfEnergy ();
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy);
  void fillPdf1d();

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  double QValue[999];				// ENDF parameter and Q values from file 3
  double sigDiff;					// precision for cross-section reconstruction
  std::vector<int> MtNumbers;				// MT numbers temporary
  rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  std::vector<int> nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  std::vector<int> nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  std::vector<int> nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  rowd energyFile5;
  rowd energyPdfFile5;
  rowd energyCdfFile5;
  rowd ein, cdf, pdf;
  matrixd2 cdf2d, pdf2d,ein2d;
  matrixd3 cdf3d, pdf3d;
#ifdef USE_ROOT
  TRandom *fRnd;
#endif
  ClassDef(TNudyEndfEnergy, 1) // class for an ENDF energy reconstruction
};
#endif
