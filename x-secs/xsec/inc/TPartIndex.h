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
#include <TDatabasePDG.h>

#define DICLEN 12        // Number of process cross sections 
#define FNPROC 17        // Number of total processes
#define FNPART 464       // Total number of particles

class TPartIndex: public TObject {

public:
   static TPartIndex *I() {if(!fgPartIndex) fgPartIndex=new TPartIndex();
      return fgPartIndex;}
   TPartIndex();
   virtual ~TPartIndex();

   // Process index <- G4 process*1000+subprocess
   Int_t ProcIndex(Short_t proccode) const;      
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
   void SetPartTable(const Int_t *PDG, Int_t np);
   
   // PDG code <- G5 particle number
   Int_t PDG(Int_t i) const {return fPDG[i];}
   // PDG code <- particle name 
   Int_t PDG(const Char_t* pname) const;
   // Particle name <- G5 particle number
   const Char_t *PartName(Int_t i) const 
   {return TDatabasePDG::Instance()->GetParticle(fPDG[i])->GetName();}
   // G5 particle index <- PDG code
   Int_t PartIndex(Int_t pdg) const {Int_t np=fNPart; 
      while(np--) if(fPDG[np]==pdg) break; return np;}
   // G5 particle index <- particle name 
   Int_t PartIndex(const Char_t *partname) const {
      return PartIndex(PDG(partname));}
   // Number of particles 
   Int_t NPart() const {return fNPart;}
 
    // Number of particles with reactions
   void SetNPartReac(Int_t np) {fNpReac=np;}
   Int_t NPartReac() const {return fNpReac;}
   TDatabasePDG *DBPdg() const {return fDBPdg;}

   void Print(Option_t *option="") const;

private:
   static TPartIndex *fgPartIndex;

   static const Int_t   fNProc=FNPROC;    // Number of processes
   static const char   *fPrName[FNPROC];  // Process name
   static const Short_t fPCode[FNPROC];   // G4 process codes

   Int_t    fNPart;         // Total number of particles
   Int_t   *fPDG;           // [fNPart] PDG code of all part

   Int_t fNpReac;           // Number of particles with reactions

   TDatabasePDG *fDBPdg;    // Pointer to the augmented pdg database

   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
