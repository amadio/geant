#ifndef TNudyEndfNuPh_H
#define TNudyEndfNuPh_H

#include "TNudyEndfRecoPoint.h"

#define PI acos(-1.0)

class  TNudyEndfNuPh : public TNudyEndfRecoPoint {

public: 
  TNudyEndfNuPh ();
  TNudyEndfNuPh (TNudyEndfFile *file);
  virtual ~TNudyEndfNuPh ();
private:
  double recursionLinearNuPh(double x1, double x2, double sig1, double sig2, std::vector<double> x, std::vector<double> sig);
  double A, AWR, ABN, QX;                         // standard ENDF parameters
  int NR, NP;                         // standard ENDF parameters for range and interpolation
  std::vector<double> eintFile1,nutFile1,einFile1,nuFile1;
  std::vector<double> eindFile1,nudFile1,einphFile1,phFile1;
  std::vector<double> einfFile1,heatFile1;
  std::vector<double> cnc;
  std::vector<int> nbt1,int1;
#ifdef USE_ROOT
  ClassDef(TNudyEndfNuPh, 1) // class for an ENDF reconstruction
#endif
};
#endif
