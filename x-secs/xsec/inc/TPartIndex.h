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

#define DICLEN 12

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPartIndex                                                           //
//                                                                      //
// Particle index singleton for various particle translation functions  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TObject.h>

#define FNPROC 17
#define FNPREA 53
#define NMAT 118
#define FNPART 464

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
   Int_t NPart() const {return fNPart;}
   Int_t Rcode(const Char_t* reac) const {Int_t nr=fNProc; 
      while(nr--) if(!strcmp(reac,fPName[nr])) break; if(nr<0) return nr;
		     return fPCode[nr];}
   Int_t NReac() const {return fNProc;}
   const Char_t *ReacName(Int_t i) const {return fPName[i];}
   Int_t Rindex(const Char_t *reac) const {Int_t nr=fNProc;
      while(nr--) if(!strcmp(reac,fPName[nr])) break; return nr;}
   static const char* MatSymb(Int_t i) {return fMatSymbol[i-1];}
   static const char* MatName(Int_t i) {return fMatName[i-1];}
   static Float_t WMat(Int_t z) {return fWmate[z];}
   void Print(Option_t *option="") const;

private:
   static TPartIndex *fgPartIndex;

   static const Int_t   fNmat=NMAT;       // Number of Materials
   static const Char_t *fMatSymbol[NMAT]; // Symbol of Material
   static const Char_t *fMatName[NMAT];   // Name of Material
   static const Float_t fWmate[NMAT];     // Weight of a mole in grams

   static const Int_t   fNProc=FNPROC;    // Number of processes
   static const char   *fPName[FNPROC];   // Process name
   static const Short_t fPCode[FNPROC];   // G4 process codes

   Int_t    fNPart;         // Total number of particles
   Int_t    fNPart30;       // length of the char store
   Int_t   *fPDG;           // [fNPart] PDG code of all part
   Char_t  *fPnames;        // [fNPart30] Names of all part

   Int_t    fNReac;          // Number of particles with reactions
   Short_t  fPDGReac[FNPREA];// Correspondence PDG <-> particle number with reac

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
