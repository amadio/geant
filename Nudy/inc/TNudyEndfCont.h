#ifndef ROOT_TNudyEndfCont
#define ROOT_TNudyEndfCont

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
//class TNudyEndfSec;
#include "TNudyEndfRecord.h"

class TNudyEndfCont: public TNudyEndfRecord {
public:
  TNudyEndfCont();
  TNudyEndfCont(double c1, double c2,
		Int_t l1, Int_t l2, Int_t n1, Int_t n2);

  virtual void SetCont(double c1, double c2,
		Int_t l1, Int_t l2, Int_t n1, Int_t n2);

  virtual double GetC1() const {return fC1;}
  virtual double GetC2() const {return fC2;}
  virtual Int_t    GetL1() const {return fL1;}
  virtual Int_t    GetL2() const {return fL2;}
  virtual Int_t    GetN1() const {return fN1;}
  virtual Int_t    GetN2() const {return fN2;}

  void DumpENDF(Int_t mat, Int_t mf, Int_t mt, Int_t& ns,Int_t flags);
  static Char_t * F2F(double f, char s[]);
 protected:
  double fC1;         // C1 of the CONT record
  double fC2;         // C2 of the CONT record
  Int_t    fL1;         // L1 of the CONT record
  Int_t    fL2;         // L2 of the CONT record
  Int_t    fN1;         // N1 of the CONT record
  Int_t    fN2;         // N2 of the CONT record

  ClassDef(TNudyEndfCont,1)

};

#endif
