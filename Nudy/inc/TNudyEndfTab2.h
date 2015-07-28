#ifndef ROOT_TNudyEndfTab2
#define ROOT_TNudyEndfTab2

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
#include "TNudyEndfCont.h"
#include <RConfig.h>

class TNudyEndfTab2: public TNudyEndfCont {
public:
  TNudyEndfTab2();
  TNudyEndfTab2(Double_t c1, Double_t c2,
		   Int_t l1, Int_t l2, Int_t n1, Int_t n2);
  virtual ~TNudyEndfTab2();
  virtual void SetCont(Double_t c1, Double_t c2,
		       Int_t l1, Int_t l2, Int_t n1, Int_t n2);

  Int_t GetNR() const {return fN1;}
  Int_t GetNZ() const {return fN2;}

  Int_t GetNBT(Int_t i) const {return fNBT[i];}
  void SetNBT(Int_t iel, Int_t i) {fNBT[i]=iel;}
  Int_t GetINT(Int_t i) const {return fINT[i];}
  void SetINT(Int_t iel, Int_t i) {fINT[i]=iel;}
  void DumpENDF(Int_t mat, Int_t mf,Int_t mt,Int_t& ns,Int_t flags);
  Int_t *NBT(){return fNBT;}
  Int_t *INT(){return fINT;}

private:
  Int_t    *fNBT;       //[fN1]
  Int_t    *fINT;       //[fN1]

  ClassDef(TNudyEndfTab2,1)

};

#endif
