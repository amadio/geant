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

#include "TNudyEndfCont.h"

class TNudyEndfTab2 : public TNudyEndfCont {
public:
  TNudyEndfTab2();
  TNudyEndfTab2(double c1, double c2, int l1, int l2, int n1, int n2);
  virtual ~TNudyEndfTab2();
  virtual void SetCont(double c1, double c2, int l1, int l2, int n1, int n2);

  int GetNR() const { return fN1; }
  int GetNZ() const { return fN2; }

  int GetNBT(int i) const { return fNBT[i]; }
  void SetNBT(int iel, int i) { fNBT[i] = iel; }
  int GetINT(int i) const { return fINT[i]; }
  void SetINT(int iel, int i) { fINT[i] = iel; }
  void DumpENDF(int mat, int mf, int mt, int &ns, int flags);
  int *NBT() { return fNBT; }
  int *INT() { return fINT; }

private:
  int *fNBT; //[fN1]
  int *fINT; //[fN1]
#ifdef USE_ROOT
  ClassDef(TNudyEndfTab2, 1)
#endif
};

#endif
