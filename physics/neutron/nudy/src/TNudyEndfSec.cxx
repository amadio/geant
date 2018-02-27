/*
   This is the main class supporting an ENDF section in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

#include <TNudyEndfSec.h>
#include "Geant/TNudyEndfCont.h"
#ifdef USE_ROOT
ClassImp(TNudyEndfSec)
#endif

    //_______________________________________________________________________________
    TNudyEndfSec::TNudyEndfSec()
    : fMAT(0), fMF(0), fMT(0), fC1(0), fC2(0), fL1(0), fL2(0), fN1(0), fN2(0)
{
  //
  // Default constructor
  //
  strcpy(fName, "");
  fRecs = new TList();
}

//_______________________________________________________________________________
TNudyEndfSec::TNudyEndfSec(int mat, int mf, int mt, double c1, double c2, int l1, int l2, int n1, int n2)
    : fMAT(mat), fMF(mf), fMT(mt), fC1(c1), fC2(c2), fL1(l1), fL2(l2), fN1(n1), fN2(n2)
{
  //
  // Standard constructor
  //
  snprintf(fName, 12, "%04d-%02d-%03d", fMAT, fMF, fMT);
  fRecs = new TList();
}

//______________________________________________________________________________
TNudyEndfSec::~TNudyEndfSec()
{
  // printf("Deleting Record\n");
  if (fRecs) {
    fRecs->Delete();
    SafeDelete(fRecs);
  }
}

//______________________________________________________________________________
TNudyEndfRecord *TNudyEndfSec::GetRecord(int RecNo)
{
  //
  // RecNo is the index of the record in the section
  // RecNo = 0 is the first record after the HEAD
  // the maximum value of RecNo is the index of the last record before SEND
  //
  if (RecNo >= 0 && RecNo <= fRecs->LastIndex())
    return (TNudyEndfRecord *)fRecs->At(RecNo);
  else {
    Error("TNudyEndfSec::GetRecord(int)", "Could not find record %d on tape", RecNo);
    return NULL;
  }
}

//______________________________________________________________________________
void TNudyEndfSec::DumpENDF(int flags = 1)
{
  // HEAD
  char s1[14], s2[14];
  TNudyEndfCont::F2F(fC1, s1);
  TNudyEndfCont::F2F(fC2, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fL1, fL2, fN1, fN2);
  printf("%4d%2d%3d%5d", fMAT, fMF, fMT, 1);
  if (flags)
    printf("  ---HEAD\n");
  else
    printf("\n");
  // RECORDS
  int ns = 2;
  for (int i = 0; i <= fRecs->LastIndex(); i++)
    ((TNudyEndfRecord *)fRecs->At(i))->DumpENDF(fMAT, fMF, fMT, ns, flags);
  // SEND
  printf("%11s%11s%11d%11d%11d%11d%4d%2d%3d%5d", " 0.000000+0", " 0.000000+0", 0, 0, 0, 0, fMAT, fMF, 0, 99999);
  if (flags)
    printf("  ---SEND\n");
  else
    printf("\n");
}
