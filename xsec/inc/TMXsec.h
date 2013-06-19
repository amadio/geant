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

#include <TNamed.h>
#include <TEXsec.h>

class TMXsec : public TNamed {

public:
   TMXsec();
   TMXsec(const Char_t* name, const Char_t *title,
	  const Int_t z[], const Int_t a[], const Float_t w[], 
	  Int_t nel, Float_t dens, Bool_t weight=kFALSE);
   ~TMXsec() {}
   Float_t Xlength(Int_t part, Float_t en);
   TEXsec *SampleInt(Int_t part, Double_t en, Int_t &reac);

private:
   Int_t          fNEbins;        // number of energy bins
   Double_t       fEmin;          // min tab energy
   Double_t       fEmax;          // max tab energy
   Double_t       fEDelta;        // multiplicative energy delta
   Double_t       fElDelta;       // logarithmic energy delta

   Int_t         fNElems;  // Number of elements
   TEXsec      **fElems; // [fNElems] List of elements composing this material
   Float_t      *fTotXL;   // Total x-sec for this material
   Float_t      *fRelXS;   // Relative x-sec for this material

   ClassDef(TMXsec,1)  //Material X-secs

};


#endif
