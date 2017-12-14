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

#include "TNudyEndfRecord.h"

class TNudyEndfCont : public TNudyEndfRecord {
public:
  TNudyEndfCont();
  TNudyEndfCont(double c1, double c2, int l1, int l2, int n1, int n2);

  virtual void SetCont(double c1, double c2, int l1, int l2, int n1, int n2);

  virtual double GetC1() const { return fC1; }
  virtual double GetC2() const { return fC2; }
  virtual int GetL1() const { return fL1; }
  virtual int GetL2() const { return fL2; }
  virtual int GetN1() const { return fN1; }
  virtual int GetN2() const { return fN2; }

  void DumpENDF(int mat, int mf, int mt, int &ns, int flags);
  static char *F2F(double f, char s[]);

protected:
  double fC1; // C1 of the CONT record
  double fC2; // C2 of the CONT record
  int fL1;    // L1 of the CONT record
  int fL2;    // L2 of the CONT record
  int fN1;    // N1 of the CONT record
  int fN2;    // N2 of the CONT record

  ClassDef(TNudyEndfCont, 1)
};

#endif
