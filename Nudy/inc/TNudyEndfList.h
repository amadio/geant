#ifndef ROOT_TNudyEndfList
#define ROOT_TNudyEndfList

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

class TNudyEndfList: public TNudyEndfCont {
public:
  TNudyEndfList();
  TNudyEndfList(Double_t c1, Double_t c2,
		   Int_t l1, Int_t l2, Int_t n1, Int_t n2);
  virtual ~TNudyEndfList();
  virtual void SetCont(Double_t c1, Double_t c2,
		       Int_t l1, Int_t l2, Int_t n1, Int_t n2);

  Int_t GetNPL() const {return fN1;}
  Double_t GetLIST(Int_t i) const {return fList[i];}
  void SetLIST(Double_t el, Int_t i) {fList[i]=el;}
  void DumpENDF(Int_t mat, Int_t mf, Int_t mt, Int_t& ns, Int_t flags);
private:
  Double_t *fList;      //[fN1]

  ClassDef(TNudyEndfList,1)

};

#endif
