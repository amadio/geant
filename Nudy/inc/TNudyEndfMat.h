#ifndef ROOT_TNudyEndfMat
#define ROOT_TNudyEndfMat

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF material in R-ENDF format

   Units used
   -------------------------------------------------------------------
   Quantity           Units
   energies           electron-volts (eV)
   angles             dimensionless cosines of the angle
   cross sections     barns
   temperatures       Kelvin
   mass               units of the neutron mass
   angular distrib    probability per unit-cosine
   energy distrib     probability per eV
   energy-angle dist  probability per unit-cosine per eV
   half life          seconds
   -------------------------------------------------------------------

*/

class TString;
class TNudyEndfFile;

#include <Riostream.h>
#include <TObject.h>
#include <TList.h>
#include "TNudyEndfFile.h"
#include <RConfig.h>

class TNudyEndfMat : public TObject {
public:
  TNudyEndfMat();
  TNudyEndfMat(int mat, int za, double awr, int lrp, bool lfi, int nlib, int nmod);

  virtual ~TNudyEndfMat();
  void SetName(char *name) {
    strncpy(fName, name, 11);
    fName[11] = '\0';
  }
  void SetMAT(int mat) { fMAT = mat; }
  void SetZA(int za) { fZA = za; }
  void SetAWR(double awr) { fAWR = awr; }
  void SetLRP(int lrp) { fLRP = lrp; }

  void SetLFI(bool lfi) { fLFI = lfi; }
  void SetNLIB(int nlib) { fNLIB = nlib; }
  void SetNMOD(int nmod) { fNMOD = nmod; }
  void SetELIS(double elis) { fELIS = elis; }
  void SetSTA(bool sta) { fSTA = sta; }
  void SetLIS(int lis) { fLIS = lis; }
  void SetLISO(int liso) { fLISO = liso; }
  void SetNFOR(int nfor) { fNFOR = nfor; }
  void SetAWI(double awi) { fAWI = awi; }
  void SetEMAX(double emax) { fEMAX = emax; }
  void SetLREL(int lrel) { fLREL = lrel; }
  void SetNSUB(int nsub) { fNSUB = nsub; }
  void SetNVER(int nver) { fNVER = nver; }
  void SetTEMP(double temp) { fTEMP = temp; }
  void SetLDRV(int ldrv) { fLDRV = ldrv; }
  void SetNWD(int nwd) { fNWDm5 = nwd - 5; }
  // void SetDesc   (const char *desc, int i);
  void SetDesc(const TString desc, int i);
  void SetNXC(int nxc);
  void SetZSYMAM(const char *zsymam) {
    strncpy(fZSYMAM, zsymam, 11);
    fZSYMAM[11] = '\0';
  }
  void SetALAB(const char *alab) {
    strncpy(fALAB, alab, 11);
    fALAB[11] = '\0';
  }
  void SetEDATE(const char *edate) {
    strncpy(fEDATE, edate, 10);
    fEDATE[10] = '\0';
  }
  void SetAUTH(const char *auth) {
    strncpy(fAUTH, auth, 33);
    fAUTH[33] = '\0';
  }
  void SetREF(const char *ref) {
    strncpy(fREF, ref, 21);
    fREF[21] = '\0';
  }
  void SetDDATE(const char *ddate) {
    strncpy(fDDATE, ddate, 10);
    fDDATE[10] = '\0';
  }
  void SetRDATE(const char *rdate) {
    strncpy(fRDATE, rdate, 10);
    fRDATE[10] = '\0';
  }
  void SetENDATE(int endate) { fENDATE = endate; }
  void SetHSUB(char *hsub, int i) { strncpy(fHSUB[i], hsub, 66), fHSUB[i][66] = '\0'; }
  void SetMFn(int mf, int i) { fMFn[i] = mf; }
  void SetMTn(int mt, int i) { fMTn[i] = mt; }
  void SetNCn(int nc, int i) { fNCn[i] = nc; }
  void SetMODn(int mod, int i) { fMODn[i] = mod; }

  void Add(TNudyEndfFile *file) { fFiles->Add(file); }

  const char *GetName() const { return fName; }
  int GetMAT() const { return fMAT; }
  int GetZA() const { return fZA; }
  double GetAWR() const { return fAWR; }
  int GetLRP() const { return fLRP; }
  bool GetLFI() const { return fLFI; }
  int GetNLIB() const { return fNLIB; }
  int GetNMOD() const { return fNMOD; }
  double GetELIS() const { return fELIS; }
  bool GetSTA() const { return fSTA; }
  int GetLIS() const { return fLIS; }
  int GetLISO() const { return fLISO; }
  int GetNFOR() const { return fNFOR; }
  double GetAWI() const { return fAWI; }
  double GetEMAX() const { return fEMAX; }
  int GetLREL() const { return fLREL; }
  int GetNSUB() const { return fNSUB; }
  int GetNVER() const { return fNVER; }
  double GetTEMP() const { return fTEMP; }
  int GetLDRV() const { return fLDRV; }
  int GetNWD() const { return fNWDm5 + 5; }
  const char *GetDesc(int i) const;
  int GetNXC() const { return fNXC; }
  const char *GetZSYMAM() const { return fZSYMAM; }
  const char *GetALAB() const { return fALAB; }
  const char *GetEDATE() const { return fEDATE; }
  const char *GetAUTH() const { return fAUTH; }
  const char *GetREF() const { return fREF; }
  const char *GetDDATE() const { return fDDATE; }
  const char *GetRDATE() const { return fRDATE; }
  int GetENDATE() const { return fENDATE; }
  const char *GetHSUB(int i) const { return fHSUB[i]; }
  const int *GetMFa() const { return fMFn; }
  const int *GetMTa() const { return fMTn; }
  const int *GetNCa() const { return fNCn; }
  const int *GetMODa() const { return fMODn; }
  int GetMFn(int i) const { return fMFn ? fMFn[i] : -1; }
  int GetMTn(int i) const { return fMTn ? fMTn[i] : -1; }
  int GetNCn(int i) const { return fNCn ? fNCn[i] : -1; }
  int GetMODn(int i) const { return fMODn ? fMODn[i] : -1; }

  // TList    fFiles;         //! List of the files of this material

  void Print(const char *) const;
  void DumpENDF(int flags);
  TNudyEndfFile *GetFile(int MF);
  TList *GetFiles() { return fFiles; }

private:
  char fName[12]; // Material name
  int fMAT;       // MAT number
  int fZA;        // Standard identifier ZA = 1000.0 × Z + A
  double fAWR;    // Material mass in atomic units
  int fLRP;       // True if resonance parameters given in File 2
  bool fLFI;      // True if this material is fissile
  int fNLIB;      // Library identifier
  int fNMOD;      // Modification number for this material
  double fELIS;   // Excitation energy of the target nucleus relative to 0.0 for the ground state.
  bool fSTA;      // True if target is stable
  int fLIS;       // State number of the target nucleus; 0=ground state
  int fLISO;      // Isomeric state number; 0=ground state
  int fNFOR;      // Library format.
  double fAWI;    // Mass of the projectile in neutron mass units
  double fEMAX;   // Upper limit of the energy range for evaluation
  int fLREL;      // Library release number; for example, LREL=2 for the ENDF/B-VI.2 library.
  int fNSUB;      // Sub-library number
  int fNVER;      // Library version number; for example, NVER=7 for version ENDF/B-VII.
  double fTEMP;   // Target temperature (K) for data Doppler broadened. 0 for primary evalua- tions.
  int fLDRV;      // Derived material evaluation flag
  int fNWDm5;     // Number of records with descriptive text minus 5 standard lines
  // char   (*fDesc)[66];// [fNWD-5] records with descriptive text
  TString *fDesc;    // [fNWDm5] records with descriptive text
  int fNXC;          // Number of records in the directory for this material
  char fZSYMAM[12];  // Character representation of the material  in the form Z-cc-AM
  char fALAB[12];    // Mnemonic for the originating laboratory(s)
  char fEDATE[11];   // Date of evaluation given in the form "11EVAL-DEC74"
  char fAUTH[34];    // Author(s) name(s)
  char fREF[22];     // Primary reference for the evaluation
  char fDDATE[11];   // Original distribution date given in the form "DIST-DEC74"
  char fRDATE[11];   // Date and number of the last revision to this evaluation "REV2-DEC74"
                     // where "2" in column 4 is the release number and must be equal to LREL
  int fENDATE;       // Master File entry date in the form yyyymmdd
  char fHSUB[3][67]; // Identifier for the library contained on three successive records
  int *fMFn;         // [fNXC] ENDF file number of the nth section.
  int *fMTn;         // [fNXC] ENDF reaction designation of the nth section.
  int *fNCn;         // [fNXC] Number of records in the nth section excluding the SEND record.
  int *fMODn;        // [fNXC] Modification indicator for the nth section.

  TList *fFiles; // List of the files of this material

  ClassDef(TNudyEndfMat, 1)
};

#endif
