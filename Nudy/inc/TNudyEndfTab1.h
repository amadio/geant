#ifndef ROOT_TNudyEndfTab1
#define ROOT_TNudyEndfTab1

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
class TNudyEndfTab1: public TNudyEndfCont {
public:
  TNudyEndfTab1();
  TNudyEndfTab1(TNudyEndfTab1* tab, Int_t n1, Int_t n2);
  TNudyEndfTab1(Double_t c1, Double_t c2,
		Int_t l1, Int_t l2, Int_t n1, Int_t n2);
  virtual ~TNudyEndfTab1();
  virtual void SetCont(Double_t c1, Double_t c2,
		       Int_t l1, Int_t l2, Int_t n1, Int_t n2);

  Int_t GetNR() const {return fN1;}
  Int_t GetNP() const {return fN2;}

  Int_t GetNBT(Int_t i) const {return fNBT[i];}
  void SetNBT(Int_t iel, Int_t i) {fNBT[i]=iel;}

  Int_t GetINT(Int_t i) const {return fINT[i];}
  void SetINT(Int_t iel, Int_t i) {fINT[i]=iel;}

  Double_t GetX(Int_t i) const {return fX[i];}
  void SetX(Double_t x, Int_t i) {fX[i]=x;}

  void Equate(TNudyEndfTab1 *tab);

  Double_t* X(){return fX;}
  Double_t* Y(){return fY;}
  Int_t* NBT(){return fNBT;}
  Int_t* INT(){return fINT;}

  Double_t GetY(Int_t i) const {return fY[i];}
  void SetY(Double_t y, Int_t i) {fY[i]=y;}

  void DumpENDF(Int_t mat, Int_t mf,Int_t mt,Int_t& ns,Int_t flags);

private:
  Int_t    *fNBT;       //[fN1]
  Int_t    *fINT;       //[fN1]
  Double_t *fX;         //[fN2]
  Double_t *fY;         //[fN2]

  ClassDef(TNudyEndfTab1,1)

};

#endif
