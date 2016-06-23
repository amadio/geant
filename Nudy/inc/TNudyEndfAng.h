#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include "TNudyEndfRecoPoint.h"

class  TNudyEndfAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfAng ();
  TNudyEndfAng (TNudyEndfFile *file);
  virtual ~TNudyEndfAng ();
private:
  double recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  std::vector<double>ein,cdf,pdf,lCoef1;
  std::vector<std::vector<double> >pdf2d,cdf2d,lCoef;
  std::vector<double> cosFile4;
  std::vector<double> cosPdfFile4;
  std::vector<double> cosCdfFile4;
  std::vector<int> nbt1,int1;
  std::vector<int> MtNumbers;				// MT numbers
  std::vector<int> MtLct;				// LCT numbers
  int nr, np;                         // standard ENDF parameters  
#ifdef USE_ROOT
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
