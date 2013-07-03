// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TEXsec
#define ROOT_TEXsec


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEXSec                                                               //
//                                                                      //
// X-section for G5 per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
class TFile;
class TGraph;
class TList;
class TPXsec;
class TPartIndex;

#define NELEM 118         // Total number of materials

class TEXsec : public TNamed {

public:
   TEXsec();
   TEXsec(Int_t z, Int_t a, Float_t dens, Int_t np);
   ~TEXsec();
   Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nxsec);
   Bool_t AddPartXS(Int_t kpart, const Float_t xsec[], const Int_t dict[]);
   Bool_t AddPartIon(Int_t kpart, const Float_t dedx[]);
   Bool_t AddPartMS(Int_t kpart, const Float_t angle[], const Float_t ansig[],
		    const Float_t length[], const Float_t lensig[]);
   
   static const char* EleSymb(Int_t z) {return fEleSymbol[z-1];}
   static const char* EleName(Int_t z) {return fEleName[z-1];}
   static Float_t WEle(Int_t z) {return fWElem[z-1];}
   static Int_t NElem() {return fNElem;}

   Int_t Ele() const {return fEle;}
   Double_t Dens() const {return fDens;}
   Double_t Emin() const {return fEmin;}
   Double_t Emax() const {return fEmax;}
   Int_t NEbins() const {return fNEbins;}
   Double_t EilDelta() const {return fEilDelta;}
   Float_t XS(Int_t pindex, Int_t rindex, Float_t en) const;
   Float_t DEdx(Int_t pindex, Float_t en) const;
   Bool_t MS(Int_t index, Float_t en, Float_t &ang, Float_t &asig, 
		  Float_t &len, Float_t &lsig) const;
   TGraph *XSGraph(const char* part, const char *reac, 
		   Float_t emin, Float_t emax, Int_t nbin) const;
   TGraph *DEdxGraph(const char* part, 
		   Float_t emin, Float_t emax, Int_t nbin) const;
   TGraph *MSGraph(const char* part, const char *what,
		   Float_t emin, Float_t emax, Int_t nbin) const;

   Float_t Lambda(Int_t pindex, Double_t en) const;
   Int_t SampleReac(Int_t pindex, Double_t en) const;
   
   static Bool_t FloatDiff(Double_t a, Double_t b, Double_t prec) {
      return TMath::Abs(a-b)>0.5*TMath::Abs(a+b)*prec;
   }

   const Double_t *Cuts() const {return fCuts;}

   void DumpPointers() const;
   void Draw(Option_t *option);
   Bool_t Resample();

   Bool_t Prune();

   static Int_t NLdElems() {return fNLdElems;}
   static TEXsec *Element(Int_t i) {
      if(i<0 || i>=fNLdElems) return 0; return fElements[i];}

   static TEXsec *GetElement(Int_t z, Int_t a=0, TFile *f=0);

private:
   static const Int_t   fNElem=NELEM;       // Number of Elements
   static const Char_t *fEleSymbol[NELEM]; // Symbol of Element
   static const Char_t *fEleName[NELEM];   // Name of Element
   static const Float_t fWElem[NELEM];     // Weight of a mole in grams

   Int_t          fEle;     // Element code Z*10000+A*10+metastable level
   Float_t        fDens;    // Density in g/cm3
   Double_t       fAtcm3;   // Atoms per cubic cm unit density
   Double_t       fEmin;    // Minimum of the energy Grid
   Double_t       fEmax;    // Maximum of the energy Grid
   Int_t          fNEbins;  // Number of log steps in energy
   Double_t       fEilDelta; // Inverse log energy step
   const Double_t *fEGrid;  //! Common energy grid
   Int_t          fNRpart;  // Number of particles with reaction
   TPXsec        *fPXsec;   // [fNRpart] Cross section table per particle
   Double_t      *fCuts;    // [fNRpart] Just a placeholder for the moment

   static Int_t   fNLdElems; //! number of loaded elements
   static TEXsec *fElements[NELEM]; //! databases of elements

   TPartIndex *   fPartIndex; // TPartIndex class

   ClassDef(TEXsec,1)  // Element X-secs

};


#endif
