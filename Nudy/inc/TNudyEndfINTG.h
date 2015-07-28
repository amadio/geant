#ifndef ROOT_TNudyEndfINTG
#define ROOT_TNudyEndfINTG

#include <Riostream.h>
#include "TNudyEndfCont.h"

class TNudyEndfINTG: public TNudyEndfCont {

 public:

  TNudyEndfINTG();
  TNudyEndfINTG(Int_t nrow, Int_t ndigit); 
  Int_t* GetKIJ() {return fKIJ;}
  Int_t GetNdigit(){return fNdigit;}
  void SetIJ(Int_t ij[2]){fII = ij[0]; fJJ = ij[1];}
  void SetNrow(Int_t nrow) { fNrow = nrow;}
  void SetNdigit(Int_t ndigit){ fNdigit = ndigit;}
  void SetKIJ(Int_t kij[18]);
  void DumpENDF(Int_t mat, Int_t mf, Int_t mt, Int_t& ns, Int_t flags);


 private:

  Int_t fKIJ[18]; 
  Int_t fNrow;
  Int_t fNdigit;
  Int_t fII;
  Int_t fJJ;
  ClassDef(TNudyEndfINTG,1)

};

#endif
