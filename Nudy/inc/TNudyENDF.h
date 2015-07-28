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
#include <sstream>
#include <string>
#include <cstring>

/*

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
*/
/* FORTRAN Reading Routines */
/*
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
  void FGetMTF(DEFCHARD str, int mtf[4] DEFCHARL);
  void FGetCONT(DEFCHARD str, double c[2], int nl[4], int mtf[4] DEFCHARL);
  void FGetFloat(DEFCHARD str, double arr[6] DEFCHARL);
  void FGetInt(DEFCHARD str, int iarr[6] DEFCHARL);
  void FGetINTG(DEFCHARD str, int &ndigit,int ij[2], int kij[18], int mtf[4] DEFCHARL);
}

*/
#define LINLEN 256

class TNudyENDF: public TObject {
 public:
  TNudyENDF();
  TNudyENDF(const char *nFileENDF, const char *nFileRENDF, const char *opt="new",unsigned char loglev=0);
  virtual ~TNudyENDF();
   void SetLogLev(unsigned char loglev) {fLogLev=loglev;}
   unsigned char GetLogLev() const {return fLogLev;}
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
   void GetSEND(const int pmtf[3]);
   void GetFEND(const int pmtf[3]);
   void GetMEND(const int pmtf[3]);
   void GetTEND();
   void ToEndSec();
   TNudyEndfTape* GetTape(){return this->fTape;}

   void CheckSEND(const int pmtf[3]) const;
   void CheckFEND(const int pmtf[3]) const;
   void CheckMEND(const int pmtf[3]) const;
   void CheckTEND() const;

   void GetMTF(int mtf[4]) const {
     //FGetMTF(PASSCHARD(fLine), mtf PASSCHARL(fLine));
     std::string s0(fLine);
     std::string s2 = s0.substr(66,4);
     std::istringstream ss;
     ss.str(s2);
     mtf[0]=0;
     ss >> mtf[0];
     ss.str("");
     ss.clear(); 
     //std::cout << "DEBUG:: mtf[0] " << mtf[0] << "  ";
     
     s2 = s0.substr(70,2);
     ss.str(s2);
     mtf[1]=0;
     ss >> mtf[1];
     ss.str("");
     ss.clear(); 
     //std::cout << "DEBUG:: mtf[1] " << mtf[1] << "  ";

     s2 = s0.substr(72,3);
     ss.str(s2);
     mtf[2]=0;
     ss >> mtf[2];
     ss.str("");
     ss.clear(); 
     //std::cout << "DEBUG:: mtf[2] " << mtf[2] << "  ";

     s2 = s0.substr(75,5);
     ss.str(s2);
     mtf[3]=0;
     ss >> mtf[3];
     ss.str("");
     ss.clear(); 
     //std::cout << "DEBUG:: mtf[3] " << mtf[3] << ".\n";

   }

   void GetCONT(double c[2], int nl[4], int mtf[4]) const {
     //FGetCONT(PASSCHARD(fLine), c, nl, mtf PASSCHARL(fLine));
     int ii;
     //std::cout << " Process :: GetCONT::-> " << fLine << std::endl;
     std::string tmp;
     std::istringstream ss;
     std::string s0(fLine);
     std::vector<std::string> strNum(6);

     for (ii = 0; ii < 6; ii++){
       strNum[ii] = s0.substr(ii*11,11);
     }
     ss.str(""); ss.clear();


     for (ii = 0; ii <2; ii++) {
       tmp.swap(strNum[ii]);
       std::size_t alien=tmp.find_last_of("+-");
       if( 0<alien && alien!=std::string::npos) tmp.replace(alien,1,std::string("E") + tmp[alien]);
       ss.str(tmp);
       c[ii]=0.0;
       ss >> c[ii];
       ss.str(""); ss.clear(); tmp="";
     }

     for (ii = 0; ii < 4; ii++) {
       nl[ii]=0;
       //std::cout << strNum[ii+2] << std::endl;
       ss.str(strNum[ii+2]);
       ss >> nl[ii];
       ss.str(""); ss.clear(); tmp=""; 
     }
     GetMTF(mtf);
     
   }

   void GetFloat(double c[6]) const {
     //FGetFloat(PASSCHARD(fLine), c PASSCHARL(fLine));
     int ii;
     std::string tmp;
     std::vector<std::string>strNum(6);
     std::istringstream ss;
     std::string s0(fLine);
     std::string s1 = s0.substr(0,66);
     ss.str(s1);
     
     for (ii = 0; ii < 6; ii++) ss >> strNum[ii]; ss.str(""); ss.clear();

     for (ii = 0; ii < 6; ii++) {
       c[ii]=0.0;
       tmp.swap(strNum[ii]);
       std::size_t alien=tmp.find_last_of("+-");
       if( 0<alien && alien!=std::string::npos) tmp.replace(alien,1,std::string("E") + tmp[alien]);
       ss.str(tmp);
       ss >> c[ii];
       ss.str(""); ss.clear(); tmp="";
     }
   }

   void GetInt(int n[6]) const {
     //FGetInt(PASSCHARD(fLine), n PASSCHARL(fLine));
     std::string s0(fLine);
     std::string s1 = s0.substr(0,66);
     std::istringstream ss(s1);
     for (int ii = 0; ii < 6; ii++) {
       n[ii]=0;
       ss >> n[ii];
     }
   }

   void GetINTG(int ndigit, int ij[2] ,int kij[18], int mtf[4] )const{
     //FGetINTG(PASSCHARD(fLine), ndigit, ij, kij, mtf PASSCHARL(fLine));
     int ii;
     std::string s0(fLine);
     std::string s1 = s0.substr(0,66);
     std::istringstream ss(s1);

     ij[0]=0; ij[1]=0;
     ss >> ij[0]; ss >> ij[1];
     for (ii = 0; ii < 18; ii++) {
       kij[ii]=0;
       ss >> kij[ii];
     }
     GetMTF(mtf);

     /*
     switch(ndigit) {
     case 2:
       ij[0]=0; ij[1]=0;
       ss >> ij[0];     ss >> ij[1];
       for (ii = 0; ii < 18; ii++) {
	 kij[ii]=0;
	 ss >> kij[ii];
       }
       GetMTF(mtf);
     case 3:
       ij[0]=0; ij[1]=0;
       ss >> ij[0];     ss >> ij[1];
       for (ii = 0; ii < 13; ii++) {
	 kij[ii]=0;
	 ss >> kij[ii];
       }
       GetMTF(mtf);
     case 4:
       ij[0]=0; ij[1]=0;
       ss >> ij[0];     ss >> ij[1];
       for (ii = 0; ii < 11; ii++) {
	 kij[ii]=0;
	 ss >> kij[ii];
       }
       GetMTF(mtf);
     case 5:
       ij[0]=0; ij[1]=0;
       ss >> ij[0];     ss >> ij[1];
       for (ii = 0; ii < 9; ii++) {
	 kij[ii]=0;
	 ss >> kij[ii];
       }
       GetMTF(mtf);
     case 6:
       ij[0]=0; ij[1]=0;
       ss >> ij[0];     ss >> ij[1];
       for (ii = 0; ii < 8; ii++) {
	 kij[ii]=0;
	 ss >> kij[ii];
       }
       GetMTF(mtf);
     default: std::cout << " ERROR : INVALID NDIGIT\n";
     }
     */
             
   }
     
   void DumpENDF(int flags);
   
private:
   static const char fkElNam[119][4];
   static const char fkElIso[4][2];
   
   unsigned char        fLogLev;        //  Log Level Flag
   std::ifstream       fENDF;          //! Input fENDF tape
   TFile         *fRENDF;         //! Output fRENDF file
   char         fLine[LINLEN];  //! Buffer to read the line
   TNudyEndfTape *fTape;          //! Support link for the tape structure
   TNudyEndfMat  *fMat;           //! Support link for the current material
   
   ClassDef(TNudyENDF, 1) // class for an ENDF data file
};

#endif

