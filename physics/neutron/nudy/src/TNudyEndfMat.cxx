/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

*/

/*
#include <string.h>

#include <TString.h>
*/

#include <TError.h>
#include "Geant/TNudyEndfCont.h"
#include "Geant/TNudyEndfMat.h"

using namespace Nudy;

#ifdef USE_ROOT
ClassImp(TNudyEndfMat)
#endif

    //_______________________________________________________________________________
    TNudyEndfMat::TNudyEndfMat()
    : fMAT(0), fZA(0), fAWR(0), fLRP(-1), fLFI(false), fNLIB(0), fNMOD(0), fELIS(0), fSTA(true), fLIS(0), fLISO(0),
      fNFOR(0), fAWI(0), fEMAX(0), fLREL(0), fNSUB(0), fNVER(0), fTEMP(0), fLDRV(0), fNWDm5(0), fDesc(NULL), fNXC(0),
      fENDATE(0), fMFn(NULL), fMTn(NULL), fNCn(NULL), fMODn(NULL)
{
  //
  // Default constructor
  //
  strcpy(fName, "");
  strcpy(fZSYMAM, "");
  strcpy(fALAB, "");
  strcpy(fEDATE, "");
  strcpy(fAUTH, "");
  strcpy(fREF, "");
  strcpy(fDDATE, "");
  strcpy(fRDATE, "");
  fFiles = new TList();
}

//_______________________________________________________________________________
TNudyEndfMat::TNudyEndfMat(int mat, int za, double awr, int lrp, bool lfi, int nlib, int nmod)
    : fMAT(mat), fZA(za), fAWR(awr), fLRP(lrp), fLFI(lfi), fNLIB(nlib), fNMOD(nmod), fELIS(0), fSTA(true), fLIS(0),
      fLISO(0), fNFOR(0), fAWI(0), fEMAX(0), fLREL(0), fNSUB(0), fNVER(0), fTEMP(0), fLDRV(0), fNWDm5(0), fDesc(NULL),
      fNXC(0), fENDATE(0), fMFn(NULL), fMTn(NULL), fNCn(NULL), fMODn(NULL)
{
  //
  // Standard constructor
  //
  strcpy(fName, "");
  strcpy(fZSYMAM, "");
  strcpy(fALAB, "");
  strcpy(fEDATE, "");
  strcpy(fAUTH, "");
  strcpy(fREF, "");
  strcpy(fDDATE, "");
  strcpy(fRDATE, "");
  fFiles = new TList();
}

//______________________________________________________________________________
TNudyEndfMat::~TNudyEndfMat()
{
  //  printf("Destroying Mat %s\n",fName);
  if (fDesc) {
    delete[] fDesc;
    fDesc = 0;
  }
  if (fMFn) {
    delete[] fMFn;
    fMFn = 0;
  }
  if (fMTn) {
    delete[] fMTn;
    fMTn = 0;
  }
  if (fNCn) {
    delete[] fNCn;
    fNCn = 0;
  }
  if (fMODn) {
    delete[] fMODn;
    fMODn = 0;
  }
  fFiles->Delete();
  SafeDelete(fFiles);
}

//_______________________________________________________________________________
const char *TNudyEndfMat::GetDesc(int i) const
{
  if (i < 0 || i >= fNWDm5) {
    Error("GetDesc", "index %d out of bounds [%d,%d]", i, 0, fNWDm5);
    return 0;
  } else
    return fDesc[i].Data();
}

//_______________________________________________________________________________
// void TNudyEndfMat::SetDesc(const char *desc, int i)
void TNudyEndfMat::SetDesc(const TString desc, int i)
{
  if (i < 0 || i >= fNWDm5) {
    Error("GetDesc", "index %d out of bounds [%d,%d]", i, 0, fNWDm5);
  } else {
    if (!fDesc) fDesc = new TString[fNWDm5];
    fDesc[i]          = desc;
  }
}

//_______________________________________________________________________________
void TNudyEndfMat::SetNXC(int nxc)
{
  fNXC = nxc;
  if (fMFn) delete[] fMFn;
  if (fMTn) delete[] fMTn;
  if (fNCn) delete[] fNCn;
  if (fMODn) delete[] fMODn;

  fMFn  = new int[fNXC];
  fMTn  = new int[fNXC];
  fNCn  = new int[fNXC];
  fMODn = new int[fNXC];
}

//_______________________________________________________________________________
void TNudyEndfMat::DumpENDF(int flags = 1)
{
  int ns = 1;
  int i  = 0;
  // Dump what was read into this classi(file 1 mt 451)
  char s1[14], s2[14];
  TNudyEndfCont::F2F(fZA, s1);
  TNudyEndfCont::F2F(fAWR, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fLRP, fLFI, fNLIB, fNMOD);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---HEAD\n");
  else
    printf("\n");
  TNudyEndfCont::F2F(fELIS, s1);
  TNudyEndfCont::F2F(fSTA, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fLIS, fLISO, 0, fNFOR);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---CONT\n");
  else
    printf("\n");
  TNudyEndfCont::F2F(fAWI, s1);
  TNudyEndfCont::F2F(fEMAX, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fLREL, 0, fNSUB, fNVER);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---CONT\n");
  else
    printf("\n");
  TNudyEndfCont::F2F(fTEMP, s1);
  TNudyEndfCont::F2F(0.0, s2);
  printf("%11s%11s%11d%11d%11d%11d", s1, s2, fLDRV, 0, GetNWD(), fNXC);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---CONT\n");
  else
    printf("\n");
  printf("%11s%11s%10s %33s", fZSYMAM, fALAB, fEDATE, fAUTH);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---TEXT\n");
  else
    printf("\n");
  printf(" %21s%10s %10s%12s%-8d   ", fREF, fDDATE, fRDATE, " ", fENDATE);
  printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
  if (flags)
    printf("  ---TEXT\n");
  else
    printf("\n");

  for (i = 0; i < 3; i++) {
    printf("%66s", fHSUB[i]);
    printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
    if (flags)
      printf("  ---TEXT\n");
    else
      printf("\n");
  }
  // Check no of lines here
  for (i = 0; i < fNWDm5; i++) {
    printf("%66s", GetDesc(i));
    printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
    if (flags)
      printf("  ---TEXT\n");
    else
      printf("\n");
  }
  for (i = 0; i < fNXC; i++) {
    printf("%11s%11s%11d%11d%11d%11d", " ", " ", fMFn[i], fMTn[i], fNCn[i], fMODn[i]);
    printf("%4d%2d%3d%5d", fMAT, 1, 451, ns++);
    if (flags)
      printf("  ---CONT\n");
    else
      printf("\n");
  }
  printf("%11s%11s%11d%11d%11d%11d%4d%2d%3d%5d", " 0.000000+0", " 0.000000+0", 0, 0, 0, 0, fMAT, 1, 0, 99999);
  if (flags)
    printf("  ---SEND\n");
  else
    printf("\n");
  // Files
  for (i = 0; i <= fFiles->LastIndex(); i++) {
    TNudyEndfFile *file = (TNudyEndfFile *)fFiles->At(i);
    file->DumpENDF(flags);
  }
  // MEND
  printf("%66s", " 0.000000+0 0.000000+0          0          0          0          0");
  printf("%4d%2d%3d%5d", 0, 0, 0, 0);
  if (flags)
    printf("  ---MEND\n");
  else
    printf("\n");
}

//_______________________________________________________________________________
TNudyEndfFile *TNudyEndfMat::GetFile(int MF)
{
  for (int i = 0; i <= this->fFiles->LastIndex(); i++) {
    TNudyEndfFile *thisFile = (TNudyEndfFile *)this->fFiles->At(i);
    if (thisFile->GetMF() == MF) return thisFile;
  }
  Error("TNudyEndfMat::GetFile(int)", "Could not find file %d on tape", MF);
  return NULL;
}

//_______________________________________________________________________________
void TNudyEndfMat::Print(const char *op) const
{
  TString sop = op;
  if (sop.Contains("h") || sop.Contains("a")) {
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::endl << "Header for Material " << fName << "MAT " << fMAT << std::endl;
    std::cout << std::setfill('-') << std::setw(80) << "-" << std::endl;
    if (sop.Contains("F")) {
      std::cout << std::setprecision(4) << std::scientific << std::setfill(' ') << std::setw(11) << double(fZA)
                << std::setw(11) << fAWR << std::setw(11) << int(fLRP) << std::setw(11) << int(!fLFI) << std::setw(11)
                << fNLIB << std::setw(11) << fNMOD << std::endl;

      std::cout << std::setprecision(4) << std::scientific << std::setfill(' ') << std::setw(11) << fELIS
                << std::setw(11) << double(fSTA == 0) << std::setw(11) << int(fLIS) << std::setw(11) << int(fLISO)
                << std::setw(11) << 0 << std::setw(11) << fNFOR << std::endl;

      std::cout << std::setprecision(4) << std::scientific << std::setfill(' ') << std::setw(11) << fAWI
                << std::setw(11) << fEMAX << std::setw(11) << fLREL << std::setw(11) << 0 << std::setw(11) << fNSUB
                << std::setw(11) << fNVER << std::endl;

      std::cout << std::setprecision(4) << std::scientific << std::setfill(' ') << std::setw(11) << fTEMP
                << std::setw(11) << 0.0 << std::setw(11) << fLDRV << std::setw(11) << 0 << std::setw(11) << fNWDm5 + 5
                << std::setw(11) << fNXC << std::endl;

      std::cout << fZSYMAM << fALAB << fEDATE << " " << fAUTH << std::endl;

      std::cout << " " << fREF << fDDATE << " " << fRDATE << std::setw(20) << fENDATE << std::endl;

      std::cout << "----" << fHSUB[0] << std::endl
                << "-----" << fHSUB[1] << std::endl
                << "------" << fHSUB[2] << std::endl;

      for (int i = 0; i < fNWDm5; ++i)
        std::cout << fDesc[i] << std::endl;

      for (int i = 0; i < fNXC; ++i)
        std::cout << std::setw(33) << fMFn[i] << std::setw(11) << fMTn[i] << std::setw(11) << fNCn[i] << std::setw(11)
                  << fMODn[i] << std::endl;
    }
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
  }
}
