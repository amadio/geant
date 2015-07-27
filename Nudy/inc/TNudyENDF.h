#ifndef ROOT_TNudyENDF
#define ROOT_TNudyENDF

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class to read a file in ENDF format and write a file
   in R-ENDF format

*/


class TFile;
class TNudyEndfTape;
class TNudyEndfMat;
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfList;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfINTG;

#include <TObject.h>
#include <TList.h>
#include <Riostream.h>
#include <RConfig.h>

#ifndef WIN32
# define type_of_call
# define DEFCHARD     const char* 
# define DEFCHARL   , const int 
# define PASSCHARD(string) string 
# define PASSCHARL(string) , strlen(string) 
#else
# define type_of_call  _stdcall
# define DEFCHARD   const char* , const int        
# define DEFCHARL          
# define PASSCHARD(string) string, strlen(string) 
# define PASSCHARL(string) 
#endif

/* FORTRAN Reading Routines */

#ifndef WIN32
#define FGetMTF   fgetmtf_
#define FGetCONT  fgetcont_
#define FGetFloat fgetfloat_
#define FGetInt   fgetint_
#define FGetINTG  fgetintg_
#else
#define FGetMTF   FGETMTF
#define FGetCONT  FGETCONT
#define FGetFloat FGETFLOAT
#define FGetInt   FGETINT
#define FGetINTG  FGETINTG
#endif

extern "C" type_of_call {
  void FGetMTF(DEFCHARD str, Int_t mtf[4] DEFCHARL);
  void FGetCONT(DEFCHARD str, double c[2], Int_t nl[4], Int_t mtf[4] DEFCHARL);
  void FGetFloat(DEFCHARD str, double arr[6] DEFCHARL);
  void FGetInt(DEFCHARD str, Int_t iarr[6] DEFCHARL);
  void FGetINTG(DEFCHARD str, Int_t &ndigit,Int_t ij[2], Int_t kij[18], Int_t mtf[4] DEFCHARL);
}

#define LINLEN 256

class TNudyENDF: public TObject {
 public:
  TNudyENDF();
  TNudyENDF(const Char_t *nFileENDF, const Char_t *nFileRENDF, Option_t *opt="new",UChar_t loglev=0);
  virtual ~TNudyENDF();
   void SetLogLev(UChar_t loglev) {fLogLev=loglev;}
   UChar_t GetLogLev() const {return fLogLev;}
   void Process();
   void Process(TNudyEndfMat *mat);
   void Process(TNudyEndfFile *file);
   void Process(TNudyEndfSec* sec);
   void Process(TNudyEndfCont *secCont);
   void Process(TNudyEndfList *secList);
   void Process(TNudyEndfTab1 *secTab1);
   void Process(TNudyEndfTab2 *secTab2);
   void Process(TNudyEndfINTG *secINTG);
   void ProcessF1(TNudyEndfSec* sec);
   void ProcessF2(TNudyEndfSec* sec);
   void ProcessF3(TNudyEndfSec* sec);
   void ProcessF4(TNudyEndfSec* sec);
   void ProcessF5(TNudyEndfSec* sec);
   void ProcessF6(TNudyEndfSec* sec);
   void ProcessF7(TNudyEndfSec* sec);
   void ProcessF8(TNudyEndfSec* sec);
   void ProcessF9(TNudyEndfSec* sec);
   void ProcessF10(TNudyEndfSec* sec);
   void ProcessF12(TNudyEndfSec* sec);
   void ProcessF13(TNudyEndfSec* sec);
   void ProcessF14(TNudyEndfSec* sec);
   void ProcessF15(TNudyEndfSec* sec);
   void ProcessF23(TNudyEndfSec* sec);
   void ProcessF26(TNudyEndfSec* sec);
   void ProcessF27(TNudyEndfSec* sec);
   void ProcessF28(TNudyEndfSec* sec);
   void ProcessF30(TNudyEndfSec* sec);
   void ProcessF31(TNudyEndfSec* sec);
   void ProcessF32(TNudyEndfSec* sec);
   void ProcessF33(TNudyEndfSec* sec);
   void ProcessF34(TNudyEndfSec* sec);
   void ProcessF35(TNudyEndfSec* sec);
   void ProcessF40(TNudyEndfSec* sec);
   void GetSEND(const Int_t pmtf[3]);
   void GetFEND(const Int_t pmtf[3]);
   void GetMEND(const Int_t pmtf[3]);
   void GetTEND();
   void ToEndSec();
   TNudyEndfTape* GetTape(){return this->fTape;}

   void CheckSEND(const Int_t pmtf[3]) const;
   void CheckFEND(const Int_t pmtf[3]) const;
   void CheckMEND(const Int_t pmtf[3]) const;
   void CheckTEND() const;

   void GetMTF(Int_t mtf[4]) const {
     FGetMTF(PASSCHARD(fLine), mtf PASSCHARL(fLine));
   }

   void GetCONT(double c[2], Int_t nl[4], Int_t mtf[4]) const {
     FGetCONT(PASSCHARD(fLine), c, nl, mtf PASSCHARL(fLine));
   }

   void GetFloat(double c[6]) const {
     FGetFloat(PASSCHARD(fLine), c PASSCHARL(fLine));
   }

   void GetInt(Int_t n[6]) const {
     FGetInt(PASSCHARD(fLine), n PASSCHARL(fLine));
   }
   void GetINTG(Int_t ndigit, Int_t ij[2] ,Int_t kij[18], Int_t mtf[4] )const{
     FGetINTG(PASSCHARD(fLine), ndigit, ij, kij, mtf PASSCHARL(fLine));
   }
     
   void DumpENDF(Int_t flags);
   
private:
   static const char fkElNam[119][4];
   static const char fkElIso[4][2];
   
   UChar_t        fLogLev;        //  Log Level Flag
   std::ifstream       fENDF;          //! Input fENDF tape
   TFile         *fRENDF;         //! Output fRENDF file
   Char_t         fLine[LINLEN];  //! Buffer to read the line
   TNudyEndfTape *fTape;          //! Support link for the tape structure
   TNudyEndfMat  *fMat;           //! Support link for the current material
   
   ClassDef(TNudyENDF, 1) // class for an ENDF data file
};

#endif

