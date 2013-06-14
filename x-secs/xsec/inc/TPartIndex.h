// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPartIndex
#define ROOT_TPartIndex


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPartIndex                                                           //
//                                                                      //
// Particle index singleton for various particle translation functions  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>

#define DICLEN 12        // Number of process cross sections 
#define FNPROC 17        // Number of total processes
#define FNPREA 53        // Number of particles with processes
#define NMAT 118         // Total number of materials
#define FNPART 464       // Total number of particles

class TPartIndex: public TObject {

public:
   static TPartIndex *I() {if(!fgPartIndex) fgPartIndex=new TPartIndex();
      return fgPartIndex;}
   TPartIndex();
   virtual ~TPartIndex();

   Short_t ProcIndex(Short_t proccode) const;
   const Char_t* ProcNameCode(Int_t proccode) const;
   const Char_t* ProcNameIndex(Int_t procindex) const;
   Int_t ProcCode(Int_t procindex) const {return fPCode[procindex];}
   Short_t NProc() const {return fNProc;}

   void SetPartTable(char **names, Int_t *PDG, Int_t np);
   
   Int_t PDG(Int_t i) const {return fPDG[i];}
   Int_t PDG(const Char_t* pname) const {Int_t nr=fNPart;
      while(nr--) if(!strcmp(pname,&fPnames[30*nr])) break;
      if(nr<0) return -12345678; return fPDG[nr];}
   const Char_t *PartName(Int_t i) const {return &fPnames[30*i];}

   Int_t PartIndex(Int_t pdg) const {Int_t np=fNPart; 
      while(np--) if(fPDG[np]==pdg) break; return np;}
   Int_t PartIndex(const Char_t *partname) const {
      return PartIndex(PDG(partname));}

   Int_t NPart() const {return fNPart;}

   Int_t NPartReac() const {return fNpReac;}
   Bool_t SetPartReac(const Int_t reacpdg[], Int_t np) {
      if(np != fNpReac) {
	 Error("SetPartReac","np is %d and should be %d\n",np,fNpReac);
	 exit(1);
	 return kFALSE;}
      for(Int_t i=0; i<np; ++i) fPDGReac[i]=reacpdg[i];
      return kTRUE;}
   Int_t PartReacIndex(Int_t pdg) {Int_t np=fNpReac;
      while(np--) if(fPDGReac[np]==pdg) break; return np;}

   Int_t Rcode(const Char_t* reac) const {Int_t nr=fNProc; 
      while(nr--) if(!strcmp(reac,fPrName[nr])) break; if(nr<0) return nr;
		     return fPCode[nr];}
   const Char_t *ReacName(Int_t i) const {return fPrName[i];}
   Int_t Rindex(const Char_t *reac) const {Int_t nr=fNProc;
      while(nr--) if(!strcmp(reac,fPrName[nr])) break; return nr;}
   static const char* EleSymb(Int_t i) {return fEleSymbol[i-1];}
   static const char* EleName(Int_t i) {return fEleName[i-1];}
   static Float_t WEle(Int_t z) {return fWElem[z];}
   void Print(Option_t *option="") const;

private:
   static TPartIndex *fgPartIndex;

   static const Int_t   fNmat=NMAT;       // Number of Elements
   static const Char_t *fEleSymbol[NMAT]; // Symbol of Element
   static const Char_t *fEleName[NMAT];   // Name of Element
   static const Float_t fWElem[NMAT];     // Weight of a mole in grams

   static const Int_t   fNProc=FNPROC;    // Number of processes
   static const char   *fPrName[FNPROC];  // Process name
   static const Short_t fPCode[FNPROC];   // G4 process codes

   Int_t    fNPart;         // Total number of particles
   Int_t    fNPart30;       // length of the char store
   Int_t   *fPDG;           // [fNPart] PDG code of all part
   Char_t  *fPnames;        // [fNPart30] Names of all part

   Int_t    fNpReac;         // Number of particles with reactions
   Int_t    fPDGReac[FNPREA];// Correspondence PDG <-> particle number with reac

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
