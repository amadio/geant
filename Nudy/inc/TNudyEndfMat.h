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

class TNudyEndfMat: public TObject {
public:
  TNudyEndfMat();
  TNudyEndfMat(Int_t mat, Int_t za, double awr, Char_t lrp, Bool_t lfi, Int_t nlib, Int_t nmod);
  virtual ~TNudyEndfMat();
  void SetName   (Char_t *name)  {strncpy(fName,name,11); fName[11]='\0';}
  void SetMAT    (Int_t mat)     {fMAT = mat;}
  void SetZA     (Int_t za)      {fZA = za;}
  void SetAWR    (double awr)  {fAWR = awr;}
  void SetLRP    (Char_t lrp)    {fLRP = lrp;}
  void SetLFI    (Bool_t lfi)    {fLFI = lfi;}
  void SetNLIB   (Int_t nlib)    {fNLIB = nlib;}
  void SetNMOD   (Int_t nmod)    {fNMOD = nmod;}
  void SetELIS   (double elis) {fELIS = elis;}
  void SetSTA    (Bool_t sta)    {fSTA = sta;}
  void SetLIS    (UChar_t lis)   {fLIS = lis;}
  void SetLISO   (UChar_t liso)  {fLISO = liso;}
  void SetNFOR   (Int_t nfor)    {fNFOR = nfor;}
  void SetAWI    (double awi)  {fAWI = awi;}
  void SetEMAX   (double emax) {fEMAX = emax;}
  void SetLREL   (Int_t lrel)    {fLREL = lrel;}
  void SetNSUB   (Int_t nsub)    {fNSUB = nsub;}
  void SetNVER   (Int_t nver)    {fNVER = nver;}
  void SetTEMP   (double temp) {fTEMP = temp;}
  void SetLDRV   (Int_t ldrv)    {fLDRV = ldrv;}
  void SetNWD    (Int_t nwd)     {fNWDm5 = nwd-5;}
  void SetDesc   (const Char_t *desc, Int_t i);
  void SetNXC    (Int_t nxc);
  void SetZSYMAM (const Char_t *zsymam)  {strncpy(fZSYMAM,zsymam,11);fZSYMAM[11]='\0';}
  void SetALAB   (const Char_t *alab)    {strncpy(fALAB,alab,11);fALAB[11]='\0';}
  void SetEDATE  (const Char_t *edate)   {strncpy(fEDATE,edate,10);fEDATE[10] = '\0';}
  void SetAUTH   (const Char_t *auth)    {strncpy(fAUTH,auth,33);fAUTH[33] = '\0';}
  void SetREF    (const Char_t *ref)     {strncpy(fREF,ref,21);fREF[21] = '\0';}
  void SetDDATE  (const Char_t *ddate)   {strncpy(fDDATE,ddate,10);fDDATE[10] = '\0';}
  void SetRDATE  (const Char_t *rdate)   {strncpy(fRDATE,rdate,10);fRDATE[10] = '\0';}
  void SetENDATE (Int_t endate)  {fENDATE = endate;}
  void SetHSUB   (Char_t *hsub, Int_t i) {strncpy(fHSUB[i],hsub,66), fHSUB[i][66] = '\0';}
  void SetMFn    (Int_t mf,  Int_t i) {fMFn[i]=mf;}
  void SetMTn    (Int_t mt,  Int_t i) {fMTn[i]=mt;}
  void SetNCn    (Int_t nc,  Int_t i) {fNCn[i]=nc;}
  void SetMODn   (Int_t mod, Int_t i) {fMODn[i]=mod;}

  void Add(TNudyEndfFile *file) {fFiles->Add(file);}
 
  const Char_t* GetName()   const {return fName;}
  Int_t         GetMAT()    const {return fMAT;}
  Int_t         GetZA()     const {return fZA;}
  double      GetAWR()    const {return fAWR;}
  Char_t        GetLRP()    const {return fLRP;}
  Bool_t        GetLFI()    const {return fLFI;}
  Int_t         GetNLIB()   const {return fNLIB;}
  Int_t         GetNMOD()   const {return fNMOD;}
  double      GetELIS()   const {return fELIS;}
  Bool_t        GetSTA()    const {return fSTA;}
  UChar_t       GetLIS()    const {return fLIS;}
  UChar_t       GetLISO()   const {return fLISO;}
  Int_t         GetNFOR()   const {return fNFOR;}
  double      GetAWI()    const {return fAWI;}
  double      GetEMAX()   const {return fEMAX;}
  Int_t         GetLREL()   const {return fLREL;}
  Int_t         GetNSUB()   const {return fNSUB;}
  Int_t         GetNVER()   const {return fNVER;}
  double      GetTEMP()   const {return fTEMP;}
  Int_t         GetLDRV()   const {return fLDRV;}
  Int_t         GetNWD()    const {return fNWDm5+5;}
  const Char_t* GetDesc(Int_t i)   const;
  Int_t         GetNXC()    const {return fNXC;}
  const Char_t* GetZSYMAM() const {return fZSYMAM;}
  const Char_t* GetALAB()   const {return fALAB;}
  const Char_t* GetEDATE()  const {return fEDATE;}
  const Char_t* GetAUTH()   const {return fAUTH;}
  const Char_t* GetREF()    const {return fREF;}
  const Char_t* GetDDATE()  const {return fDDATE;}
  const Char_t* GetRDATE()  const {return fRDATE;}
  Int_t         GetENDATE() const {return fENDATE;}
  const Char_t* GetHSUB(Int_t i) const {return fHSUB[i];}
  const Int_t*  GetMFa  () const {return fMFn;}
  const Int_t*  GetMTa  () const {return fMTn;}
  const Int_t*  GetNCa  () const {return fNCn;}
  const Int_t*  GetMODa () const {return fMODn;}
  Int_t         GetMFn  (Int_t i) const {return fMFn?fMFn[i]:-1;}
  Int_t         GetMTn  (Int_t i) const {return fMTn?fMTn[i]:-1;}
  Int_t         GetNCn  (Int_t i) const {return fNCn?fNCn[i]:-1;}
  Int_t         GetMODn (Int_t i) const {return fMODn?fMODn[i]:-1;}

  //TList    fFiles;         //! List of the files of this material

  void Print(const Option_t*) const;
  void DumpENDF(Int_t flags);
  TNudyEndfFile* GetFile(Int_t MF);
  TList* GetFiles() {return fFiles;}
private:
  Char_t   fName[12];   // Material name
  Int_t    fMAT;        // MAT number
  Int_t    fZA;         // Standard identifier ZA = 1000.0 × Z + A
  double fAWR;        // Material mass in atomic units
  Char_t   fLRP;        // True if resonance parameters given in File 2
  Bool_t   fLFI;        // True if this material is fissile
  Int_t    fNLIB;       // Library identifier 
  Int_t    fNMOD;       // Modification number for this material
  double fELIS;       // Excitation energy of the target nucleus relative to 0.0 for the ground state.
  Bool_t   fSTA;        // True if target is stable
  UChar_t  fLIS;        // State number of the target nucleus; 0=ground state
  UChar_t  fLISO;       // Isomeric state number; 0=ground state
  Int_t    fNFOR;       // Library format.
  double fAWI;        // Mass of the projectile in neutron mass units
  double fEMAX;       // Upper limit of the energy range for evaluation
  Int_t    fLREL;       // Library release number; for example, LREL=2 for the ENDF/B-VI.2 library.
  Int_t    fNSUB;       // Sub-library number
  Int_t    fNVER;       // Library version number; for example, NVER=7 for version ENDF/B-VII.
  double fTEMP;       // Target temperature (K) for data Doppler broadened. 0 for primary evalua- tions.
  Int_t    fLDRV;       // Derived material evaluation flag
  Int_t    fNWDm5;      // Number of records with descriptive text minus 5 standard lines
  //  Char_t   (*fDesc)[66];// [fNWD-5] records with descriptive text
  TString *fDesc;       // [fNWDm5] records with descriptive text
  Int_t    fNXC;        // Number of records in the directory for this material
  Char_t   fZSYMAM[12]; // Character representation of the material  in the form Z-cc-AM
  Char_t   fALAB[12];   // Mnemonic for the originating laboratory(s)
  Char_t   fEDATE[11];  // Date of evaluation given in the form "11EVAL-DEC74"
  Char_t   fAUTH[34];   // Author(s) name(s)
  Char_t   fREF[22];    // Primary reference for the evaluation
  Char_t   fDDATE[11];  // Original distribution date given in the form "DIST-DEC74"
  Char_t   fRDATE[11];  // Date and number of the last revision to this evaluation "REV2-DEC74"
                        // where "2" in column 4 is the release number and must be equal to LREL
  Int_t    fENDATE;     // Master File entry date in the form yyyymmdd
  Char_t   fHSUB[3][67];// Identifier for the library contained on three successive records
  Int_t    *fMFn;       // [fNXC] ENDF file number of the nth section.
  Int_t    *fMTn;       // [fNXC] ENDF reaction designation of the nth section.
  Int_t    *fNCn;       // [fNXC] Number of records in the nth section excluding the SEND record.
  Int_t    *fMODn;      // [fNXC] Modification indicator for the nth section.

  TList    *fFiles;      // List of the files of this material

  ClassDef(TNudyEndfMat,1)

};

#endif
