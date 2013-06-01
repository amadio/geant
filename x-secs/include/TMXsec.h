// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TMXsec
#define ROOT_TMXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>
#include <TObject.h>
class TPXsec;

class TMXsec : public TObject {

public:
   TMXsec();
   TMXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen, Int_t np);
   ~TMXsec();

private:
   Short_t        fMat; // Material code Z*10000+A*10+metastable level
   Float_t        fEmin; // Min en in GeV
   Float_t        fEmax; // Max en in Gev
   Short_t        fNen;  // Number of log steps in energy
   Double_t       fElDelta; // Log energy step
   Short_t        fNpart; // Number of particles
   TPXsec       **fPXsec; // [fNpart] Cross section table per particle
   Double_t      *fCuts; // Just a placeholder for the moment
};

#endif
