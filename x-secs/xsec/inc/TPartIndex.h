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
#define FNPART 464       // Total number of particles

class TPartIndex: public TObject {

public:
   static TPartIndex *I() {if(!fgPartIndex) fgPartIndex=new TPartIndex();
      return fgPartIndex;}
   TPartIndex();
   virtual ~TPartIndex();

   // Process index <- G4 process*1000+subprocess
   Short_t ProcIndex(Short_t proccode) const;      
   // Process name <- process index
   const Char_t* ProcNameIndex(Int_t procindex) const;
   // Process index <- G4 process*1000+subprocess
   Int_t ProcCode(Int_t procindex) const {return fPCode[procindex];}
   // Process name <- G4 process*1000+subprocess
   const Char_t* ProcNameCode(Int_t proccode) const;
   // Number of processes -- process number fNProc()-1 is the total x-sec
   // Process code <- Process name
   Int_t ProcCode(const Char_t* reac) const {Int_t nr=fNProc; 
      while(nr--) if(!strcmp(reac,fPrName[nr])) break; if(nr<0) return nr;
		     return fPCode[nr];}
   // Process index <- Process name
   Int_t ProcIndex(const Char_t *reac) const {Int_t nr=fNProc;
      while(nr--) if(!strcmp(reac,fPrName[nr])) break; return nr;}
   Short_t NProc() const {return fNProc;}
 
   // Fill the particle table
   void SetPartTable(char **names, Int_t *PDG, Int_t np);
   
   // PDG code <- G5 particle number
   Int_t PDG(Int_t i) const {return fPDG[i];}
   // PDG code <- particle name 
   Int_t PDG(const Char_t* pname) const {Int_t nr=fNPart;
      while(nr--) if(!strcmp(pname,&fPnames[30*nr])) break;
      if(nr<0) return -12345678; return fPDG[nr];}
   // Particle name <- G5 particle number
   const Char_t *PartName(Int_t i) const {return &fPnames[30*i];}
   // G5 particle index <- PDG code
   Int_t PartIndex(Int_t pdg) const {Int_t np=fNPart; 
      while(np--) if(fPDG[np]==pdg) break; return np;}
   // G5 particle index <- particle name 
   Int_t PartIndex(const Char_t *partname) const {
      return PartIndex(PDG(partname));}
   // Number of particles 
   Int_t NPart() const {return fNPart;}
 
   // Set particles with reaction
   Bool_t SetPartReac(const Int_t reacpdg[], Int_t np);

   // Number of particles with reactions
   Int_t NPartReac() const {return fNpReac;}
   Int_t PartReacIndex(Int_t index) const {return fReacPart[index];}
   const Int_t* PartReac() const {return fReacPart;}

   void Print(Option_t *option="") const;

private:
   static TPartIndex *fgPartIndex;


   static const Int_t   fNProc=FNPROC;    // Number of processes
   static const char   *fPrName[FNPROC];  // Process name
   static const Short_t fPCode[FNPROC];   // G4 process codes

   Int_t    fNPart;         // Total number of particles
   Int_t    fNPart30;       // length of the char store
   Int_t   *fPDG;           // [fNPart] PDG code of all part
   Int_t   *fReacPart;      // [fNPart] number of particel with reaction
   Char_t  *fPnames;        // [fNPart30] Names of all part

   Int_t fNpReac;           // Number of particles with reactions

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
