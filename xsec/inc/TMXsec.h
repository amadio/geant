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
class GeantTrack_v;

class TMXsec : public TNamed {

public:
   TMXsec();
   TMXsec(const Char_t* name, const Char_t *title,
	  const Int_t z[], const Int_t a[], const Float_t w[], 
	  Int_t nel, Float_t dens, Bool_t weight=kFALSE);
   ~TMXsec();
   Float_t Xlength(Int_t part, Float_t en);
   Bool_t Xlength_v(Int_t npart, const Int_t part[], const Float_t en[], Double_t lam[]);
   Float_t DEdx(Int_t part, Float_t en, Int_t &elemindx);
   Bool_t DEdx_v(Int_t npart, const Int_t part[], const Float_t en[], Float_t de[]);
   Float_t Range(Int_t part, Float_t en, Int_t &elemindx);
   Float_t Range(Int_t part, Float_t en);
   Double_t InvRange(Int_t part, Double_t step);

   void Eloss(Int_t ntracks, GeantTrack_v &tracks);
   void	ProposeStep(Int_t ntracks, GeantTrack_v &tracks, Int_t tid);
   void	SampleInt(Int_t ntracks, GeantTrack_v &tracksin, Int_t tid);
   Float_t MS(Int_t ipart, Float_t energy);

   TEXsec *SampleInt(Int_t part, Double_t en, Int_t &reac);
   Int_t SampleElement(Int_t tid); // based on # atoms/vol. for the prototype
   Int_t SampleElement(); // based on # atoms/vol. for Geant4 with tab.phys. 

   Int_t SelectElement(Int_t pindex, Int_t rindex, Double_t energy);

   static Bool_t Prune();
   void Print(Option_t * opt="") const;

private:
   TMXsec(const TMXsec&);      // Not implemented
   TMXsec& operator=(const TMXsec&);      // Not implemented

   Int_t           fNEbins;    // number of energy bins
   Int_t           fNTotXL;    // dimension of fTotXL
   Int_t           fNCharge;   // dimension of tables for charged particles
   Int_t           fNRelXS;    // dimension of fRelXS
   Double_t        fEilDelta;  // logarithmic energy delta
   const Double_t *fEGrid;     //! Energy grid

   Int_t           fNElems;    // Number of elements
   TEXsec        **fElems;     // [fNElems] List of elements composing this material
   Float_t        *fTotXL;     // [fNTotXL] Total x-sec for this material
   Float_t        *fRelXS;     // [fNRelXS] Relative x-sec for this material
   Float_t        *fDEdx;      // [fNCharge] Ionisation energy loss for this material
   Float_t        *fMSangle;   // [fNCharge] table of MS average angle
   Float_t        *fMSansig;   // [fNCharge] table of MS sigma angle
   Float_t        *fMSlength;  // [fNCharge] table of MS average lenght correction
   Float_t        *fMSlensig;  // [fNCharge] table of MS sigma lenght correction
   Double_t       *fRatios;    // [fNElems]  relative #atoms/volume; normalized
   Float_t        *fRange;     // [fNCharge] ranges of the particle in this material

   ClassDef(TMXsec,1)  //Material X-secs

};


#endif
