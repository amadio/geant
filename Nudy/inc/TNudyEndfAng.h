#ifndef TNudyEndfAng_H
#define TNudyEndfAng_H

#include "TNudyEndfRecoPoint.h"
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowint > matrixint;
typedef std::vector<rowd > matrixd2;
typedef std::vector<std::vector<rowd > > matrixd3;
typedef std::vector<std::vector<std::vector<rowd > > > matrixd4;

#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom;
#endif

class  TNudyEndfAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfAng ();
  TNudyEndfAng (TNudyEndfFile *file);
  virtual double GetCos4(int elemid, int mt, double energyK);
  virtual ~TNudyEndfAng ();
private:
  double recursionLinearLeg(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);
  void fillPdf1d();
  void fillPdf2d();
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  rowd ein,cdf,pdf,lCoef1;
  matrixd2 pdf2d,cdf2d,lCoef,ein2d;
  matrixd3 pdf3d,cdf3d;
  rowd cosFile4;
  rowd cosPdfFile4;
  rowd cosCdfFile4;
  rowint nbt1,int1;
  rowint MtNumbers;				// MT numbers
  rowint MtLct;				// LCT numbers
  int nr, np;                         // standard ENDF parameters  
#ifdef USE_ROOT
  TRandom *fRnd;
#endif
#ifdef USE_ROOT
  ClassDef(TNudyEndfAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
