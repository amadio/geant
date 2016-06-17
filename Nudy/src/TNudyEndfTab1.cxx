/*
   This is the main class supporting an ENDF section in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

#include <string.h>
#include <stdio.h>
#include <TString.h>
#include <TNudyEndfTab1.h>

//_______________________________________________________________________________
TNudyEndfTab1::TNudyEndfTab1()
  : TNudyEndfCont(), fNBT(NULL), fINT(NULL), fX(NULL), fY(NULL) {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
TNudyEndfTab1::TNudyEndfTab1(TNudyEndfTab1 *tab, int n1, int n2)
    : TNudyEndfCont(tab->GetC1(), tab->GetC2(), tab->GetL1(), tab->GetL2(), n1, n2), fNBT(new int[n1]),
      fINT(new int[n1]), fX(new double[n2]), fY(new double[n2]) {
  memcpy(fNBT, tab->NBT(), sizeof(int) * tab->GetN1());
  memcpy(fINT, tab->INT(), sizeof(int) * tab->GetN1());
  memcpy(fX, tab->X(), sizeof(double) * tab->GetN2());
  memcpy(fY, tab->Y(), sizeof(double) * tab->GetN2());
}

//_______________________________________________________________________________
TNudyEndfTab1::TNudyEndfTab1(double c1, double c2, int l1, int l2, int n1, int n2)
    : TNudyEndfCont(c1, c2, l1, l2, n1, n2), fNBT(new int[n1]), fINT(new int[n1]), fX(new double[n2]),
      fY(new double[n2]) {
  //
  // Standard constructor
  //
}

//_______________________________________________________________________________
TNudyEndfTab1::~TNudyEndfTab1() {
  // printf("Deleting Tab1\n");
  delete [] fNBT;
  delete [] fINT;
  delete [] fX;
  delete [] fY;
}

//_______________________________________________________________________________
void TNudyEndfTab1::Equate(TNudyEndfTab1 *tab) {
  SetCont(tab->GetC1(), tab->GetC2(), tab->GetL1(), tab->GetL2(), tab->GetN1(), tab->GetN2());
  memcpy(fNBT, tab->NBT(), sizeof(int) * tab->GetN1());
  memcpy(fINT, tab->INT(), sizeof(int) * tab->GetN1());
  memcpy(fX, tab->X(), sizeof(double) * tab->GetN2());
  memcpy(fY, tab->Y(), sizeof(double) * tab->GetN2());
}

//______________________________________________________________________________
void TNudyEndfTab1::SetCont(double c1, double c2, int l1, int l2, int n1, int n2) {
  TNudyEndfCont::SetCont(c1, c2, l1, l2, n1, n2);
  delete[] fNBT;
  delete[] fINT;
  delete[] fX;
  delete[] fY;
  fNBT = new int[n1];
  fINT = new int[n1];
  fX = new double[n2];
  fY = new double[n2];
}

//
// Dumps Tab1 records to ENDF
//______________________________________________________________________________
void TNudyEndfTab1::DumpENDF(int mat, int mf, int mt, int &ns, int flags) {
  char s1[14], s2[14];
  F2F(fC1, s1);
  F2F(fC2, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fL1, fL2, fN1, fN2);
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if (ns < 99999)
    ns++;
  else
    ns = 1;
  if (flags)
    printf("  ---CONT TAB1\n");
  else
    printf("\n");
  for (int i = 0; i < GetNR(); i++) { // print NBT(N) INT(N)
    if (i % 3 == 0 && i != 0) {
      printf("%4d%2d%3d%5d", mat, mf, mt, ns);
      if (ns < 99999)
        ns++;
      else
        ns = 1;
      if (flags)
        printf("  ---NBT(%d,%d,%d) TAB1\n", i - 2, i - 1, i);
      else
        printf("\n");
    }
    printf("%11d%11d", GetNBT(i), GetINT(i));
  }
  // Pad blank columns
  if (3 - (GetNR() % 3) < 3) {
    for (int i = 0; i < 3 - (GetNR() % 3); i++) {
      printf("%22s", " ");
    }
  }
  printf("%4d%2d%3d%5d", mat, mf, mt, ns);
  if (ns < 99999)
    ns++;
  else
    ns = 1;
  if (flags)
    printf("  ---NBT(%d,%d,%d) TAB1\n", (GetNR() - (GetNR() % 3)) + 1, (GetNR() - (GetNR() % 3)) + 2,
           (GetNR() - (GetNR() % 3)) + 3);
  else
    printf("\n");
  for (int i = 0; i < GetNP();) { // print 6I11
    F2F(GetX(i), s1);
    F2F(GetY(i), s2);
    printf("%11s%11s", s1, s2);
    if ((++i) % 3 == 0) {
      printf("%4d%2d%3d%5d", mat, mf, mt, ns);
      if (ns < 99999)
        ns++;
      else
        ns = 1;
      if (flags)
        printf("  ---XY(%d,%d,%d) TAB1\n", i - 2, i - 1, i);
      else
        printf("\n");
    }
  }
  // Pad Blank Columns
  if (3 - (GetNP() % 3) < 3) {
    for (int i = 0; i < 3 - (GetNP() % 3); i++) {
      printf("%22s", " ");
    }
    printf("%4d%2d%3d%5d", mat, mf, mt, ns);
    if (ns < 99999)
      ns++;
    else
      ns = 1;
    if (flags)
      printf("  ---XY(%d,%d,%d) TAB1\n", (GetNP() - (GetNP() % 3)) + 1, (GetNP() - (GetNP() % 3)) + 2,
             (GetNP() - (GetNP() % 3)) + 3);
    else
      printf("\n");
  }
}
