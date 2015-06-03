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

#include <TObject.h>
#include <TList.h>
#include <Riostream.h>
#include <RConfig.h>

class TNudyEndfTape: public TObject {
public:
   TNudyEndfTape();
   TNudyEndfTape(const Char_t *name, UChar_t loglev);
   virtual ~TNudyEndfTape();
   const Char_t* GetName() const {return fName;}

   void SetLogLev(UChar_t loglev) {fLogLev=loglev;}
   UChar_t GetLogLev() const {return fLogLev;}
   const TList* GetMats() const {return fMats;}
   void AddMat(TNudyEndfMat* mat);
   void DumpENDF(Int_t flags);
   TNudyEndfMat* GetMAT(Int_t MAT);
   TNudyEndfMat* GetMAT(Int_t Z, Int_t A);
   TNudyEndfMat* GetMAT(TString name);
private:
   UChar_t   fLogLev;    // LogLevel
   Char_t    fName[81];  // Name of the tape
   TList     *fMats;      // List of materials
  
   ClassDef(TNudyEndfTape, 1) // class for an ENDF data file
};

#endif

