#ifndef ROOT_TNudyEndfSec
#define ROOT_TNudyEndfSec

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


#include <Riostream.h>
#include <TObject.h>
#include <TList.h>
#include <RConfig.h>
#include "TNudyEndfRecord.h"



class TNudyEndfSec: public TObject {
public:
  TNudyEndfSec();
  TNudyEndfSec(Int_t mat, Int_t mf, Int_t mt, Double_t c1, Double_t c2,
		Int_t l1, Int_t l2, Int_t n1, Int_t n2);
  virtual ~TNudyEndfSec();
  const Char_t* GetName() const {return fName;}
  void Add(TNudyEndfRecord *sec) {fRecs->Add(sec);}
  TList* GetRecords(){return fRecs;}
  void DumpENDF(Int_t flags);
  TNudyEndfRecord* GetRecord(Int_t recNo);

  Int_t GetC1() const {return fC1;}
  Int_t GetC2() const {return fC2;}
  Int_t GetL1() const {return fL1;}
  Int_t GetL2() const {return fL2;}
  Int_t GetN1() const {return fN1;}
  Int_t GetN2() const {return fN2;}
  Int_t GetMAT() const {return fMAT;}
  Int_t GetMT() const {return fMT;}
  Int_t GetMF() const {return fMF;}

 private:
  Char_t   fName[12];   // Name of the section
  Short_t  fMAT;        // Mat number
  Short_t  fMF;         // File number
  Int_t    fMT;         // Section number 
  Double_t fC1;         // C1 of the HEAD record
  Double_t fC2;         // C2 of the HEAD record
  Int_t    fL1;         // L1 of the HEAD record
  Int_t    fL2;         // L2 of the HEAD record
  Int_t    fN1;         // N1 of the HEAD record
  Int_t    fN2;         // N2 of the HEAD record

  TList    *fRecs;       // List of records for this section

  ClassDef(TNudyEndfSec,1)

};

#endif
