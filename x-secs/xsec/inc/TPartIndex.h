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

#define FNPROC 16
#define FNPREA 53
#define NMAT 100

class TPartIndex: public TObject {

public:
   static TPartIndex *I() {if(!fgPartIndex) fgPartIndex=new TPartIndex();
      return fgPartIndex;}
   TPartIndex();
   virtual ~TPartIndex();
   Short_t ProcIndex(Short_t proccode) const;
   const Char_t* ProcNameCode(Int_t proccode) const;
   const Char_t* ProcNameIndex(Int_t procindex) const;
   Short_t NProc() const {return fNProc;}

   void SetPartTable(char **names, Int_t *PDG, Int_t np);
   
   Int_t PDG(Int_t i) const {return fPDG[i];}
   Int_t PDG(const Char_t* pname) const {Int_t nr=fNPart;
      while(nr--) if(!strcmp(pname,&fPnames[30*nr])) break;
      if(nr<0) return nr; return fPDG[nr];}
   const Char_t *PartName(Int_t i) const {return &fPnames[30*i];}
   Int_t PartIndex(Int_t pdg) const {Int_t np=fNPart; 
      while(np--) if(fPDG[np]==pdg) break; return np;}
   Int_t NPart() const {return fNPart;}
   Int_t Rcode(const Char_t* reac) const {Int_t nr=fNProc;
      while(nr--) if(!strcmp(reac,fPName[nr])) break; if(nr<0) return nr;
		     return fPCode[nr];}

private:
   static TPartIndex *fgPartIndex;
   static const Char_t *fMatSymbol[NMAT] //
   static const Char_t *fMatName[NMAT] // 

   Int_t    fNProc;         // Number of processes
   char    *fPName[FNPROC]; // [fNProc] Process name
   Short_t  fPCode[FNPROC]; // G4 process codes

   Int_t    fNPart;         // Total number of particles
   Int_t    fNPart30;       // length of the char store
   Short_t *fPDG;           // [fNPart] PDG code of all part
   Char_t  *fPnames;        // [fNPart30] Names of all part

   Int_t    fNReac;          // Number of particles with reactions
   Short_t  fPDGReac[FNPREA];// Correspondence PDG <-> particle number with reac

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
