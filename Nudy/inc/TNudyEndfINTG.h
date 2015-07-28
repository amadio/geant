#ifndef ROOT_TNudyEndfINTG
#define ROOT_TNudyEndfINTG

#include <Riostream.h>
#include "TNudyEndfCont.h"

class TNudyEndfINTG : public TNudyEndfCont {

public:
  TNudyEndfINTG();
  TNudyEndfINTG(int nrow, int ndigit);
  int *GetKIJ() { return fKIJ; }
  int GetNdigit() { return fNdigit; }
  void SetIJ(int ij[2]) {
    fII = ij[0];
    fJJ = ij[1];
  }
  void SetNrow(int nrow) { fNrow = nrow; }
  void SetNdigit(int ndigit) { fNdigit = ndigit; }
  void SetKIJ(int kij[18]);
  void DumpENDF(int mat, int mf, int mt, int &ns, int flags);

private:
  int fKIJ[18];
  int fNrow;
  int fNdigit;
  int fII;
  int fJJ;
  ClassDef(TNudyEndfINTG, 1)
};

#endif
