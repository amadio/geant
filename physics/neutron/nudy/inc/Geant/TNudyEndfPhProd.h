#ifndef TNudyEndfPhProd_H
#define TNudyEndfPhProd_H

#include "Geant/TNudyEndfRecoPoint.h"

namespace NudyPhysics {
class TNudyEndfRecoPoint;
}

namespace Nudy {
class TNudyEndfFile;
}

typedef std::vector<double> rowd;
typedef std::vector<int> rowi;
typedef std::vector<rowd> matrixd2;
#ifdef USE_ROOT
#include "Rtypes.h"
class TRandom3;
#endif

namespace NudyPhysics {
class TNudyEndfPhProd : public NudyPhysics::TNudyEndfRecoPoint {

public:
  TNudyEndfPhProd();
  TNudyEndfPhProd(Nudy::TNudyEndfFile *file);
  //  virtual double GetFisYield(int elemid, double energyK);
  virtual ~TNudyEndfPhProd();

private:
  int NR, NP, LP, LF;
  double ES, EG, Eint, Yk;
  rowi nbt1, int1, lpk1, lfk1, nrk1, npk1; // tab1 parameter
  rowd es, esk1, egk1, eintk1, sigPh1;     // incident energy, yield and tab1 parameters
  double A, AWR, ABN, QX;                  // standard ENDF parameters
  rowd ein, einc;                          // incident energy
#ifdef USE_ROOT
  TRandom3 *fRnd;
#endif
  ClassDef(TNudyEndfPhProd, 1) // class for an ENDF fission yield reconstruction
};

} // namespace
#endif
