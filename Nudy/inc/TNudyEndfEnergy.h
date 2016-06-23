#ifndef TNudyEndfEnergy_H
#define TNudyEndfEnergy_H

#include "TNudyEndfRecoPoint.h"

#define PI acos(-1.0)

class  TNudyEndfEnergy : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergy ();
  TNudyEndfEnergy (TNudyEndfFile *file);
  virtual ~TNudyEndfEnergy ();
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile5GenEva(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Maxwell(double x1, double x2, double pdf1, double pdf2, double energy);
  double recursionLinearFile5Watt(double x1, double x2, double pdf1, double pdf2, double energy);

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  double QValue[999];				// ENDF parameter and Q values from file 3
  double sigDiff;					// precision for cross-section reconstruction
  std::vector<int> MtNumbers;				// MT numbers temporary
  std::vector<double> fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  std::vector<int> nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  std::vector<int> nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  std::vector<int> nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  std::vector<double> energyFile5;
  std::vector<double> energyPdfFile5;
  std::vector<double> energyCdfFile5;
  std::vector<double>ein, cdf, pdf;
  std::vector<std::vector<double> >cdf2d, pdf2d;
  ClassDef(TNudyEndfEnergy, 1) // class for an ENDF energy reconstruction
};
#endif
