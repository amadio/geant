#ifndef TNudyEndfEnergyAng_H
#define TNudyEndfEnergyAng_H

#include "TNudyEndfRecoPoint.h"

#define PI acos(-1.0)
typedef std::vector<double> rowd;
typedef std::vector<int> rowint;
typedef std::vector<rowd > matrixd2;
typedef std::vector<std::vector<rowd > > matrixd3;

class  TNudyEndfEnergyAng : public TNudyEndfRecoPoint {

public: 
  TNudyEndfEnergyAng ();
  TNudyEndfEnergyAng (TNudyEndfFile *file, double []);
  virtual ~TNudyEndfEnergyAng ();
private:
  double recursionLinearFile5Prob(double x1, double x2, double pdf1, double pdf2);
  double recursionLinearFile4(int i, double x1, double x2, double pdf1, double pdf2);
  double recursionLinearProb(double x1, double x2, double pdf1, double pdf2);

  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int testv;
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  int NR2, NE;                            // standard ENDF parameters for range and interpolation
  int NR3, NE2;				// standard ENDF parameters for range and interpolation
//  double *INorm, QValue[999];				// ENDF parameter and Q values from file 3
  rowd cosFile4, energyFile5;
  rowd cosPdfFile4, energyPdfFile5;
  rowd cosCdfFile4, energyCdfFile5;
  rowd edes6,f06,r6,a6;
  rowd fE1, fP1, fE2, fP2, fE3, fP3, INorm;
  rowint MtNumbers;				// MT numbers
  rowint MtLct;				// LCT numbers
  rowint nbt1,int1;
  int nr1, np1;                         // standard ENDF parameters
  rowint nbt2,int2;
  int nr2, np2;                         // standard ENDF parameters
  rowint nbt3,int3;
  int nr3, np3;                         // standard ENDF parameters
  rowd ein,cdf,pdf,lCoef1;
  matrixd2 pdf2d,cdf2d,lCoef,ein2d;
  matrixd3 pdf3d,cdf3d;
#ifdef USE_ROOT
  ClassDef(TNudyEndfEnergyAng, 1) // class for an ENDF reconstruction
#endif
};
#endif
