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

#include "TList.h"
#include "Geant/TNudyEndfRecord.h"

namespace Nudy {
class TNudyEndfRecord;
}

namespace Nudy {

class TNudyEndfSec : public TObject {
public:
  TNudyEndfSec();
  TNudyEndfSec(int mat, int mf, int mt, double c1, double c2, int l1, int l2, int n1, int n2);
  virtual ~TNudyEndfSec();
  const char *GetName() const { return fName; }
  void Add(Nudy::TNudyEndfRecord *sec) { fRecs->Add(sec); }
  TList *GetRecords() { return fRecs; }
  void DumpENDF(int flags);
  Nudy::TNudyEndfRecord *GetRecord(int recNo);

  double GetC1() const { return fC1; }
  double GetC2() const { return fC2; }
  int GetL1() const { return fL1; }
  int GetL2() const { return fL2; }
  int GetN1() const { return fN1; }
  int GetN2() const { return fN2; }
  int GetMAT() const { return fMAT; }
  int GetMT() const { return fMT; }
  int GetMF() const { return fMF; }

private:
  char fName[12]; // Name of the section
  short fMAT;     // Mat number
  short fMF;      // File number
  int fMT;        // Section number
  double fC1;     // C1 of the HEAD record
  double fC2;     // C2 of the HEAD record
  int fL1;        // L1 of the HEAD record
  int fL2;        // L2 of the HEAD record
  int fN1;        // N1 of the HEAD record
  int fN2;        // N2 of the HEAD record

  TList *fRecs; // List of records for this section
#ifdef USE_ROOT
  ClassDef(TNudyEndfSec, 1)
#endif
};

} // namespace
#endif
