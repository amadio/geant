/*
   This is the main class supporting an ENDF section in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

#include <string.h>

#include <TString.h>
#include <TNudyEndfList.h>

ClassImp(TNudyEndfList)

    //_______________________________________________________________________________
    TNudyEndfList::TNudyEndfList()
    : TNudyEndfCont(), fList(NULL) {
  //
  // Default constructor
  //
}

//_______________________________________________________________________________
TNudyEndfList::TNudyEndfList(double c1, double c2, int l1, int l2, int n1, int n2)
    : TNudyEndfCont(c1, c2, l1, l2, n1, n2), fList(new double[n1]) {
  //
  // Standard constructor
  //
}

//______________________________________________________________________________
TNudyEndfList::~TNudyEndfList() { SafeDelete(fList); }
//_______________________________________________________________________________
void TNudyEndfList::SetCont(double c1, double c2, int l1, int l2, int n1, int n2) {
  TNudyEndfCont::SetCont(c1, c2, l1, l2, n1, n2);
  delete[] fList;
  fList = new double[n1];
}
void TNudyEndfList::DumpENDF(int mat, int mf, int mt, int &ns, int flags = 1) {
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
    printf("  ---CONT LIST\n");
  else
    printf("\n");

  for (int i = 0; i < GetNPL(); ++i) {

    F2F(fList[i], s1);
    printf("%11s", s1);
    // printf("%11s"," ");
    if ((i + 1) % 6 == 0) {
      printf("%4d%2d%3d%5d", mat, mf, mt, ns);
      if (ns < 99999)
        ns++;
      else
        ns = 1;
      if (flags)
        printf("  ---B(%d,%d,%d,%d,%d,%d) LIST\n", i - 5, i - 4, i - 3, i - 2, i - 1, i);
      else
        printf("\n");
    }
  }
  if (6 - (GetNPL() % 6) < 6) {
    for (int i = 0; i < 6 - (GetNPL() % 6); i++) {
      //      F2F(0.0,s1);
      // printf("%11s",s1);
      printf("%11s", " ");
    }
    printf("%4d%2d%3d%5d", mat, mf, mt, ns);
    if (ns < 99999)
      ns++;
    else
      ns = 1;
    if (flags)
      printf("  ---B(%d,%d,%d,%d,%d,%d) LIST\n", GetNPL() - (GetNPL() % 6) + 1, GetNPL() - (GetNPL() % 6) + 2,
             GetNPL() - (GetNPL() % 6) + 3, GetNPL() - (GetNPL() % 6) + 4, GetNPL() - (GetNPL() % 6) + 5,
             GetNPL() - (GetNPL() % 6) + 6);
    else
      printf("\n");
  }
}
