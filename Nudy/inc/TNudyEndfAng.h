#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include "TNudyEndfRecoPoint.h"
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd > matrixd2;

class  TNudyEndfAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfAng ();
  TNudyEndfAng (TNudyEndfFile *file);
  virtual ~TNudyEndfAng ();
private:
  double recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  void fillPdf1d();
  void fillPdf2d();
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  rowd ein,cdf,pdf,lCoef1;
  matrixd2 pdf2d,cdf2d,lCoef;
  rowd cosFile4;
  rowd cosPdfFile4;
  rowd cosCdfFile4;
  rowint nbt1,int1;
  rowint MtNumbers;				// MT numbers
  rowint MtLct;				// LCT numbers
  int nr, np;                         // standard ENDF parameters  
#ifdef USE_ROOT
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
