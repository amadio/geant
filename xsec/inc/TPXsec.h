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

#include <TObject.h>
#include <TPartIndex.h>

class TPXsec: public TObject {
public:
   TPXsec();
   TPXsec(Int_t pdg, Int_t nen, Int_t nxsec, 
	  Float_t emin, Float_t emax);
   ~TPXsec();
   void Print(Option_t *opt="") const;
   Bool_t SetPart(Int_t pdg, Int_t nen, Int_t nxsec, Float_t emin, Float_t emax);
   Bool_t SetPart(Int_t pdg, Int_t nxsec);
   Bool_t SetPartXS(const Float_t xsec[], const Short_t dict[]);
   Bool_t SetPartIon(const Float_t dedx[]);
   Bool_t SetPartMS(const Float_t angle[], const Float_t ansig[],
		    const Float_t length[], const Float_t lensig[]);
   Int_t PDG() const {return fPDG;}
   Float_t XS(Short_t rindex, Float_t en) const;
   Float_t DEdx(Float_t en) const;
   Bool_t MS(Float_t en, Float_t &ang, Float_t &asig, 
	     Float_t &len, Float_t &lsig) const;
   Int_t SampleReac(Double_t en) const;
   void Dump() const;
private:
   Int_t           fPDG;           // particle pdg code
   Int_t           fNEbins;        // number of energy bins
   Int_t           fNCbins;        // number of energy bins for dEdx and MS
   Int_t           fNXsec;         // number of reactions
   Int_t           fNTotXs;        // tot size of fXSecs
   Double_t        fEmin;          // min tab energy
   Double_t        fEmax;          // max tab energy
   Double_t        fEilDelta;      // logarithmic energy delta
   const Double_t *fEGrid;         //![fNEbins] energy grid
   Float_t        *fMSangle;       // [fNCbins] table of MS average angle
   Float_t        *fMSansig;       // [fNCbins] table of MS sigma angle
   Float_t        *fMSlength;      // [fNCbins] table of MS average lenght correction
   Float_t        *fMSlensig;      // [fNCbins] table of MS sigma lenght correction
   Float_t        *fdEdx;          // [fNCbins] table of dE/dx
   Float_t        *fTotXs;         // [fNEbins] table of total x-sec
   Float_t        *fXSecs;         // [fNTotXs] table of partial x-sec
   Short_t         fRdict[FNPROC]; // reaction dictionary from reaction number to position
                                  // in the X-sec array
   Short_t         fRmap[FNPROC];  // reaction map, from reaction position in the X-sec
                                  // array to the raction number

   ClassDef(TPXsec,1)  //Particle X-secs
};

#endif
