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

//class TString;
//class TNudyEndfSec;

#include "TNudyEndfSec.h"
#include <Riostream.h>
#include <TObject.h>
#include <TList.h>
#include <RConfig.h>

class TNudyEndfFile: public TObject {
public:
  TNudyEndfFile();
  TNudyEndfFile(Int_t mat, Int_t mf);
  virtual ~TNudyEndfFile();

  const Char_t* GetName()   const {return fName;}
  Int_t         GetMAT()    const {return fMAT;}
  Int_t         GetMF()     const {return fMF;}
  void Add(TNudyEndfSec *sec) {fSecs->Add(sec);}
  
  //TList    fMF;         //! List of the files of this material

  void DumpENDF(Int_t flags);
  TNudyEndfSec* GetSec(Int_t MT);
  TList* GetSections(){return fSecs;}
private:
  Char_t   fName[9];    // File Name
  Int_t    fMAT;        // MAT number
  Int_t    fMF;         // MF number

  TList    *fSecs;         // List of the sections of this file

  ClassDef(TNudyEndfFile,1)

};

#endif
