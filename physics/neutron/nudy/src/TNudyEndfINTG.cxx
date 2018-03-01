#include "Geant/TNudyEndfINTG.h"

#ifdef USE_ROOT
ClassImp(TNudyEndfINTG)
#endif

    //______________________________________________________________________________
    TNudyEndfINTG::TNudyEndfINTG()
    : TNudyEndfCont(), fNrow(18), fNdigit(2)
{
  for (int i = 0; i < 18; fKIJ[i++] = 0)
    ;
  fII = fJJ = 0;
}

//______________________________________________________________________________
TNudyEndfINTG::TNudyEndfINTG(int nrow, int ndigit)
{
  fNrow   = nrow;
  fNdigit = ndigit;
  for (int i = 0; i < 18; fKIJ[i++] = 0)
    ;
  fII = fJJ = 0;
}

//______________________________________________________________________________
void TNudyEndfINTG::SetKIJ(int kij[18])
{
  for (int i = 0; i < 18; i++)
    fKIJ[i]  = kij[i];
}

//______________________________________________________________________________
void TNudyEndfINTG::DumpENDF(int mat, int mf, int mt, int &ns, int flags = 1)
{
  printf("%5d%5d", fII, fJJ);
  if (fNdigit >= 2 && fNdigit <= 6) printf(" ");
  int limit = 18;
  if (fNdigit == 2)
    limit = 18;
  else if (fNdigit == 3)
    limit = 13;
  else if (fNdigit == 4)
    limit = 11;
  else if (fNdigit == 5)
    limit = 9;
  else if (fNdigit == 6)
    limit = 8;
  for (int i = 0; i < limit; i++) {
    if (fKIJ[i] != 0) {
      if (fNdigit == 2)
        printf("%3d", fKIJ[i]);
      else if (fNdigit == 3)
        printf("%4d", fKIJ[i]);
      else if (fNdigit == 4)
        printf("%5d", fKIJ[i]);
      else if (fNdigit == 5)
        printf("%6d", fKIJ[i]);
      else if (fNdigit == 6)
        printf("%7d", fKIJ[i]);
    } else {
      if (fNdigit == 2)
        printf("%3s", " ");
      else if (fNdigit == 3)
        printf("%4s", " ");
      else if (fNdigit == 4)
        printf("%5s", " ");
      else if (fNdigit == 5)
        printf("%6s", " ");
      else if (fNdigit == 6)
        printf("%7s", " ");
    }
  }
  if (fNdigit == 2 || fNdigit == 5) {
    printf(" ");
  } else if (fNdigit == 3) {
    printf("   ");
  }
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if (ns < 99999)
    ns++;
  else
    ns = 1;
  if (flags)
    printf("  ---NDIGIT=%d INTG\n", fNdigit);
  else
    printf("\n");
}
