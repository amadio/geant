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

#include "TNamed.h"
#include "TPartIndex.h"
class TFile;
class TGHorizontalFrame;
class TGListBox;
class TGMainFrame;
class TGraph;
class TList;
class TPXsec;
class TRootEmbeddedCanvas;

class TEXsec : public TNamed {
public:

   enum {kCutGamma, kCutElectron, kCutPositron, kCutProton};

   TEXsec();
   TEXsec(Int_t z, Int_t a, Float_t dens, Int_t np);
   ~TEXsec();
   Bool_t AddPart(Int_t kpart, Int_t pdg, Int_t nxsec);
   Bool_t AddPartXS(Int_t kpart, const Float_t xsec[], const Int_t dict[]);
   Bool_t AddPartIon(Int_t kpart, const Float_t dedx[]);
   Bool_t AddPartMS(Int_t kpart, const Float_t angle[], const Float_t ansig[],
		    const Float_t length[], const Float_t lensig[]);
   
   Int_t Ele() const {return fEle;}
   Int_t Index() const {return fIndex;}
   void SetIndex(Int_t index) { fIndex = index; }
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
   Bool_t Lambda_v(Int_t npart, const Int_t pindex[], const Double_t en[], Double_t lam[]) const;
   Bool_t Lambda_v(Int_t npart, Int_t pindex, const Double_t en[], Double_t lam[]) const;
   Int_t SampleReac(Int_t pindex, Double_t en) const;
   Int_t SampleReac(Int_t pindex, Double_t en, Double_t randn) const;
   
   static Bool_t FloatDiff(Double_t a, Double_t b, Double_t prec) {
      return TMath::Abs(a-b)>0.5*TMath::Abs(a+b)*prec;
   }

   const Float_t *Cuts() const {return fCuts;}
   Bool_t SetCuts(const Double_t cuts[4]) {
      for(Int_t jc=0; jc<4; ++jc) fCuts[jc]=cuts[jc]; return kTRUE;}

   void DumpPointers() const;
   void Draw(Option_t *option);
   void Viewer(); // *MENU*
   void UpdateReactions();
   void SelectAll();
   void DeselectAll();
   void PreDraw();
   void ResetFrame();


   Bool_t Resample();

   Bool_t Prune();

   static Int_t NLdElems() {return fNLdElems;}
   static TEXsec *Element(Int_t i) {
      if(i<0 || i>=fNLdElems) return 0; return fElements[i];}

   static TEXsec *GetElement(Int_t z, Int_t a=0, TFile *f=0);
   static TEXsec **GetElements() {return fElements;}

private:
   TEXsec(const TEXsec &); // Not implemented
   TEXsec& operator=(const TEXsec &); // Not implemented

   Int_t          fEle;     // Element code Z*10000+A*10+metastable level
   Int_t	  fIndex;   // Index of this in TTabPhysMgr::fElemXsec 
   Float_t        fDens;    // Density in g/cm3
   Double_t       fAtcm3;   // Atoms per cubic cm unit density
   Double_t       fEmin;    // Minimum of the energy Grid
   Double_t       fEmax;    // Maximum of the energy Grid
   Int_t          fNEbins;  // Number of log steps in energy
   Double_t       fEilDelta; // Inverse log energy step
   const Double_t *fEGrid;  //! Common energy grid
   Int_t          fNRpart;  // Number of particles with reaction
   TPXsec        *fPXsec;   // [fNRpart] Cross section table per particle
   Float_t        fCuts[4]; // Production cuts "a la G4"

   static Int_t   fNLdElems; //! number of loaded elements
   static TEXsec *fElements[NELEM]; //! databases of elements

   static TGMainFrame         *fMain;           //! Main window
   static TGHorizontalFrame   *fSecond;         //! Window for the graph and the bar on left
   static TRootEmbeddedCanvas *fCanvas;         //! For the graphs
   static TGListBox           *fReactionBox;    //! Reaction list
   static TGListBox           *fParticleBox;    //! Particle list 


   ClassDef(TEXsec,2)  // Element X-secs

};


#endif
