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

#include "Geant/TNudyEndfCont.h"

class TNudyEndfTab1 : public TNudyEndfCont {
public:
  TNudyEndfTab1();
  TNudyEndfTab1(TNudyEndfTab1 *tab, int n1, int n2);
  TNudyEndfTab1(double c1, double c2, int l1, int l2, int n1, int n2);
  virtual ~TNudyEndfTab1();
  virtual void SetCont(double c1, double c2, int l1, int l2, int n1, int n2);

  int GetNR() const { return fN1; }
  int GetNP() const { return fN2; }

  int GetNBT(int i) const { return fNBT[i]; }
  void SetNBT(int iel, int i) { fNBT[i] = iel; }

  int GetINT(int i) const { return fINT[i]; }
  void SetINT(int iel, int i) { fINT[i] = iel; }

  double GetX(int i) const { return fX[i]; }
  void SetX(double x, int i) { fX[i] = x; }

  void Equate(TNudyEndfTab1 *tab);

  double *X() { return fX; }
  double *Y() { return fY; }
  int *NBT() { return fNBT; }
  int *INT() { return fINT; }

  double GetY(int i) const { return fY[i]; }
  void SetY(double y, int i) { fY[i] = y; }

  void DumpENDF(int mat, int mf, int mt, int &ns, int flags);

private:
  int *fNBT;  //[fN1]
  int *fINT;  //[fN1]
  double *fX; //[fN2]
  double *fY; //[fN2]

#ifdef USE_ROOT
  ClassDef(TNudyEndfTab1, 1)
#endif
};

#endif
