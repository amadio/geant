#ifndef TNudyEndfDoppler_H
#define TNudyEndfDoppler_H

#include "TNudyEndfRecoPoint.h"

class TNudyEndfDoppler : public TNudyEndfRecoPoint {

public:
  TNudyEndfDoppler();
  TNudyEndfDoppler(double isigDiff, double aw, double t1, double t2, std::vector<double> &x1, std::vector<double> &x2);
  virtual ~TNudyEndfDoppler();
  std::vector<double> energy;
  std::vector<double> sigma;

private:
  double OVSQPI = 0.564189583547756279;
  double boltz  = 8.617385E-05;
  double ONE = 1.0, TWO = 2.0, THH = 3.0 / TWO, TW3 = TWO / 3.0, ZERO = 0.0, HALF = ONE / TWO;
  double ZLIMI = 4;
  double tk;
  double ALPHA;
  double awri;
  double ZKT, Y2;
  double ZK2P, ZK22P, EXPAP, F0K2P, F1K2P, F2K2P, F3K2P, F4K2P;
  double ZK2, ZK22, EXPA, F0K2, F1K2, F2K2, F3K2, F4K2;
  double ZK1, F0K1, F1K1, F2K1, F3K1, F4K1;
  double FACT, AK, CK, CKY, CKY2;
  double FTAIL1 = 0, FTAIL2 = 0;
  double E1, E2, S1, S2;
  double ZK1P, F0K1P, F1K1P, F2K1P, F3K1P, F4K1P;
  double RATLOW = 1;
  double RATHIG = 2;
  double RATHLF;
  double HTEST;
  double RATMAX;
  double EMAX;
  double EMIN;
  double Y;
  double XSUM;
  int ncrs, IPP, KPP, size;
#ifdef USE_ROOT
  ClassDef(TNudyEndfDoppler, 1)
#endif
};
#endif
