// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TXsec
#define ROOT_TXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TXSec                                                                //
//                                                                      //
// X-section for G5                                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>

class TXsec {

public:
   TXsec();
   TXsec(Int_t z, Int_t a, Float_t emin, Float_t emax, Int_t nen);
   ~TXsec();
   void AddXsec(Int_t part, Int_t reac, const Float_t xsec[]) {}

private:
   struct         TXtemp {
     Short_t      fPart; // PDG code
     Short_t      fReac; // Reaction code
     Float_t     *fXs; // cross section
     TXtemp      *fNext; // Next x-section
   };
     
   Short_t        fMat; // Material code Z*1000+A
   Float_t        fEmin; // Min en in GeV
   Float_t        fEmax; // Max en in Gev
   Short_t        fNen;  // Number of log steps in energy
   Double_t       fElDelta; // Log energy step
   Float_t       *fXsec; // Cross sections
   Short_t        fNpart; // Number of particles
   Int_t         *fPcode; // Pointer to start of particle x-secs
   Short_t       *fRcode; // reaction codes for particles
   UChar_t        fNXsec; // Number of x-secs per part (+total)
   TXtemp        *fXtemp; // Temporary storage for x-secs
};

#endif
