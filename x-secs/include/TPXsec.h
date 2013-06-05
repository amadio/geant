// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPXsec
#define ROOT_TPXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPXSec                                                               //
//                                                                      //
// X-section for G5 per particle                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <RTypes.h>
#include <TObject.h>

class TPXsec: public TObject {
public:
   TPXsec();
   TPXsec(Int_t pdg, Int_t nen, Int_t nxsec, 
	  Float_t emin, Float_t emax);
   ~TPXsec();
   void Print(Option_t *opt="") const;
private:
   Int_t          fPDG;          // particle pdg code
   Int_t          fNen;          // number of energy bins
   Int_t          fNXsec;        // number of reactions
   Float_t        fEmin;         // min tab energy
   Float_t        fEmax;         // max tab energy
   Double_t       fEDelta;       // logarithmic energy delta
   Float_t       *fMSangle;      // [fNen] table of MS angle
   Float_t       *fMSlength;     // [fNen] table of MS lenght correction
   Float_t       *fdEdx;         // [fNen] table of dE/dx
   Float_t       *fTotXs;        // [fNen] table of total x-sec
   Float_t       *fXSecs;        // [fNen*fNXsec] table of partial x-sec
   Int_t         *fRdict;        // [fNXsec] reaction dictionary

   ClassDef(TPXsec,1)  //Particle X-secs
};

#endif
