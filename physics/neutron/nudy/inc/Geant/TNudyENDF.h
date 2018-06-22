#ifndef TNudyENDF_H
#define TNudyENDF_H

// Author: F.Carminati 02/05/09

/*
   This is the main class to read a file in ENDF format and write a file
   in R-ENDF format

*/

class TFile;

namespace Nudy {
class TNudyEndfTape;
class TNudyEndfMat;
class TNudyEndfFile;
class TNudyEndfSec;
class TNudyEndfCont;
class TNudyEndfList;
class TNudyEndfTab1;
class TNudyEndfTab2;
class TNudyEndfINTG;
}

#ifdef USE_ROOT
#include "Rtypes.h"
#endif

#include <sstream>
#include <fstream>
using std::ifstream;

#define LINLEN 256

namespace Nudy {

class TNudyENDF {
public:
  TNudyENDF();
  TNudyENDF(const char *nFileENDF, const char *nFileRENDF, const char *opt = "new", unsigned char loglev = 0);
  virtual ~TNudyENDF();
  bool sub = false;
  void SetEndfSub(std::string strENDFSUB)
  {
    this->sub     = true;
    this->ENDFSUB = strENDFSUB;
  }
  std::string GetEndfSubName() const { return ENDFSUB; }
  void SetLogLev(unsigned char loglev) { fLogLev = loglev; }
  void SetPreProcess(int x1) { fPrepro = x1; }
  unsigned char GetLogLev() const { return fLogLev; }
  bool GetLFI() { return fLFI; }
  void SetLFI(bool keywd) { fLFI = keywd; }
  void Process();
  void Process(Nudy::TNudyEndfMat *mat);
  void Process(Nudy::TNudyEndfFile *file);
  void Process(Nudy::TNudyEndfSec *sec);
  void Process(Nudy::TNudyEndfCont *secCont);
  void Process(Nudy::TNudyEndfList *secList);
  void Process(Nudy::TNudyEndfTab1 *secTab1);
  void Process(Nudy::TNudyEndfTab2 *secTab2);
  void Process(Nudy::TNudyEndfINTG *secINTG);
  void ProcessF1(Nudy::TNudyEndfSec *sec);
  void ProcessF2(Nudy::TNudyEndfSec *sec);
  void ProcessF3(Nudy::TNudyEndfSec *sec);
  void ProcessF4(Nudy::TNudyEndfSec *sec);
  void ProcessF5(Nudy::TNudyEndfSec *sec);
  void ProcessF6(Nudy::TNudyEndfSec *sec);
  void ProcessF7(Nudy::TNudyEndfSec *sec);
  void ProcessF8(Nudy::TNudyEndfSec *sec);
  void ProcessF9(Nudy::TNudyEndfSec *sec);
  void ProcessF10(Nudy::TNudyEndfSec *sec);
  void ProcessF12(Nudy::TNudyEndfSec *sec);
  void ProcessF13(Nudy::TNudyEndfSec *sec);
  void ProcessF14(Nudy::TNudyEndfSec *sec);
  void ProcessF15(Nudy::TNudyEndfSec *sec);
  void ProcessF23(Nudy::TNudyEndfSec *sec);
  void ProcessF26(Nudy::TNudyEndfSec *sec);
  void ProcessF27(Nudy::TNudyEndfSec *sec);
  void ProcessF28(Nudy::TNudyEndfSec *sec);
  void ProcessF30(Nudy::TNudyEndfSec *sec);
  void ProcessF31(Nudy::TNudyEndfSec *sec);
  void ProcessF32(Nudy::TNudyEndfSec *sec);
  void ProcessF33(Nudy::TNudyEndfSec *sec);
  void ProcessF34(Nudy::TNudyEndfSec *sec);
  void ProcessF35(Nudy::TNudyEndfSec *sec);
  void ProcessF40(Nudy::TNudyEndfSec *sec);
  void GetSEND(const int pmtf[3]);
  void GetFEND(const int pmtf[3]);
  void GetMEND(const int pmtf[3]);
  void GetTEND();
  void ToEndSec();
  Nudy::TNudyEndfTape *GetTape() { return this->fTape; }

  void CheckSEND(const int pmtf[3]) const;
  void CheckFEND(const int pmtf[3]) const;
  void CheckMEND(const int pmtf[3]) const;
  void CheckTEND() const;

  void GetMTF(int mtf[4]) const
  {
    std::string s0(fLine);
    std::string s2 = s0.substr(66, 4);
    std::istringstream ss;
    ss.str(s2);
    mtf[0] = 0;
    ss >> mtf[0];
    ss.str("");
    ss.clear();

    s2 = s0.substr(70, 2);
    ss.str(s2);
    mtf[1] = 0;
    ss >> mtf[1];
    ss.str("");
    ss.clear();

    s2 = s0.substr(72, 3);
    ss.str(s2);
    mtf[2] = 0;
    ss >> mtf[2];
    ss.str("");
    ss.clear();

    s2 = s0.substr(75, 5);
    ss.str(s2);
    mtf[3] = 0;
    ss >> mtf[3];
    ss.str("");
    ss.clear();
  }

  void GetCONT(double c[2], int nl[4], int mtf[4]) const
  {
    int ii;
    std::string tmp;
    std::istringstream ss;
    std::string s0(fLine);
    std::vector<std::string> strNum(6);

    for (ii = 0; ii < 6; ii++) {
      strNum[ii] = s0.substr(ii * 11, 11);
    }
    ss.str("");
    ss.clear();
    for (ii = 0; ii < 2; ii++) {
      tmp.swap(strNum[ii]);
      if (fPrepro == 0) {
        std::size_t alien = tmp.find_last_of("+-");
        if (0 < alien && alien != std::string::npos) tmp.replace(alien, 1, std::string("E") + tmp[alien]);
      }
      ss.str(tmp);
      c[ii] = 0.0;
      ss >> c[ii];
      ss.str("");
      ss.clear();
      tmp = "";
    }

    for (ii = 0; ii < 4; ii++) {
      nl[ii] = 0;
      ss.str(strNum[ii + 2]);
      ss >> nl[ii];
      ss.str("");
      ss.clear();
      tmp = "";
    }
    GetMTF(mtf);
  }

  void GetFloat(double c[6]) const
  {
    int ii;
    std::string tmp;
    std::vector<std::string> strNum(6);
    std::istringstream ss;
    std::string s0(fLine);
    std::string s1 = s0.substr(0, 66);

    for (ii      = 0; ii < 6; ii++)
      strNum[ii] = s1.substr(ii * 11, 11);

    for (ii = 0; ii < 6; ii++) {
      c[ii] = 0.0;
      tmp.swap(strNum[ii]);
      std::size_t found = tmp.find("E");
      if (fPrepro == 0 && found==std::string::npos) {
        std::size_t alien = tmp.find_last_of("+-");
        if (0 < alien && alien != std::string::npos) tmp.replace(alien, 1, std::string("E") + tmp[alien]);
      }
      ss.str(tmp);
      ss >> c[ii];
      ss.str("");
      ss.clear();
      tmp = "";
    }
  }

  void GetInt(int n[6]) const
  {
    std::string s0(fLine);
    std::string s1 = s0.substr(0, 66);
    std::istringstream ss(s1);
    for (int ii = 0; ii < 6; ii++) {
      n[ii] = 0;
      ss >> n[ii];
    }
  }

  void GetINTG(int /*ndigit*/, int ij[2], int kij[18], int mtf[4]) const
  {
    int ii;
    std::string s0(fLine);
    std::string s1 = s0.substr(0, 66);
    std::istringstream ss(s1);

    ij[0] = 0;
    ij[1] = 0;
    ss >> ij[0];
    ss >> ij[1];
    for (ii = 0; ii < 18; ii++) {
      kij[ii] = 0;
      ss >> kij[ii];
    }
    GetMTF(mtf);
  }

  void DumpENDF(int flags);

private:
  static const char fkElNam[119][4];
  static const char fkElIso[4][2];
  unsigned char fLogLev;      //  Log Level Flag
  char fLine[LINLEN];         //! Buffer to read the line
  std::string ENDFSUB;
  int fPrepro;
  bool fLFI;
  ifstream fENDF;             //! Input fENDF tape
  TFile *fRENDF;              //! Output fRENDF file
  Nudy::TNudyEndfTape *fTape; //! Support link for the tape structure
  Nudy::TNudyEndfMat *fMat;   //! Support link for the current material
#ifdef USE_ROOT
  ClassDef(TNudyENDF, 1) // class for an ENDF data file
#endif
};

} // namespace
#endif
