#ifndef ROOT_TNudyEndfTape
#define ROOT_TNudyEndfTape

// @(#)root/meta:$Id$
// Author: F.Carminati 02/05/09

/*
   This is the main class to read a file in ENDF format and write a file
   in R-ENDF format

*/

class TFile;
class TNudyEndfMat;
class TNudyENDF;

#include <TNamed.h>
#include <TList.h>
#include <Riostream.h>
#include <RConfig.h>

class TNudyEndfTape : public TNamed {
public:
  TNudyEndfTape();
  TNudyEndfTape(const char *name, unsigned char loglev);
  virtual ~TNudyEndfTape();

  void SetLogLev(unsigned char loglev) { fLogLev = loglev; }
  unsigned char GetLogLev() const { return fLogLev; }
  const TList *GetMats() const { return fMats; }
  void AddMat(TNudyEndfMat *mat);
  void DumpENDF(int flags);
  TNudyEndfMat *GetMAT(int MAT);
  TNudyEndfMat *GetMAT(int Z, int A);
  TNudyEndfMat *GetMAT(TString name);

private:
  unsigned char fLogLev; // LogLevel
  TList *fMats;    // List of materials

  ClassDef(TNudyEndfTape, 1) // class for an ENDF data file
};

#endif
