#ifndef ROOT_TNudyEndfFile
#define ROOT_TNudyEndfFile

// @(#)root/meta:$Id: TNuEndf.h 29000 2009-06-15 13:53:52Z rdm $
// Author: F.Carminati 02/05/09

/*
   This is the main class supporting an ENDF File in R-ENDF format

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

#include "Geant/TNudyEndfSec.h"
#include <TList.h>

namespace Nudy {
class TNudyEndfSec;
}

namespace Nudy {

class TNudyEndfFile : public TObject {
public:
  TNudyEndfFile();
  TNudyEndfFile(int mat, int mf);
  virtual ~TNudyEndfFile();

  const char *GetName() const { return fName; }
  int GetMAT() const { return fMAT; }
  int GetMF() const { return fMF; }
  void Add(Nudy::TNudyEndfSec *sec) { fSecs->Add(sec); }
  void SetMF(int mf);
  // TList    fMF;         //! List of the files of this material

  void DumpENDF(int flags);
  Nudy::TNudyEndfSec *GetSec(int MT);
  TList *GetSections() { return fSecs; }

private:
  char fName[9]; // File Name
  int fMAT;      // MAT number
  int fMF;       // MF number
  TList *fSecs;  // List of the sections of this file

  ClassDef(TNudyEndfFile, 1)
};

} // namespace
#endif
