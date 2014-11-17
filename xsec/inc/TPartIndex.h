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

#include "TObject.h"
#include "TDatabasePDG.h"
#include <map>

#define DICLEN 12        // Number of process cross sections 
#define FNPROC 18        // Number of total processes
#define FNPART 464       // Total number of particles
#define NELEM 118        // Total number of materials

 
enum G5proc {kTransport, kMultScatt, kIonisation, kDecay, kinElastic,
	     kElastic, kRestCapture ,kBrehms, kPairProd, kAnnihilation,
	     kCoulombScatt, kPhotoel, kCompton, kConversion, kCapture,
	     kKiller, kTotal};

class TPartIndex {

public:
   static TPartIndex *I() {if(!fgPartIndex) fgPartIndex=new TPartIndex();
      return fgPartIndex;}
   TPartIndex();
   virtual ~TPartIndex();

   // Database version
   Int_t Version() const {return fVersion;}
   Int_t VersionMajor() const {return fVersion/1000/1000;}
   Int_t VersionMinor() const {return fVersion/1000-VersionMajor()*1000;}
   Int_t VersionSub() const {return fVersion-VersionMajor()*1000000-VersionMinor()*1000;}

   // Process name <- process index
   const Char_t* ProcName(Int_t proc) const;
   // Process index <- Process name
   Int_t ProcIndex(const Char_t *reac) const {Int_t nr=fgNProc;
      while(nr--) if(!strcmp(reac,fgPrName[nr])) break; return nr;}


   // Process index <- G4 process*1000+subprocess
   Int_t ProcIndex(Int_t proccode) const;      
   // Process index <- G4 process*1000+subprocess
   static Int_t ProcCode(Int_t procindex) /* const */ {return fgPCode[procindex];}

   static Short_t NProc() /* const */ {return fgNProc;}
 
   // Fill the particle table
   void SetPartTable(const Int_t *vpdg, Int_t np);
   
   // PDG code <- G5 particle number
   Int_t PDG(Int_t i) const {return fPDG[i];}
   // PDG code <- particle name 
   Int_t PDG(const Char_t* pname) const;
   // Particle name <- G5 particle number
   const Char_t *PartName(Int_t i) const 
   {return TDatabasePDG::Instance()->GetParticle(fPDG[i])->GetName();}

   // G5 particle index <- PDG code
   Int_t PartIndex(Int_t pdg) const;

   // G5 particle index <- particle name 
   Int_t PartIndex(const Char_t *partname) const {
      return PartIndex(PDG(partname));}
   // Number of particles 
   Int_t NPart() const {return fNPart;}
 
    // Number of particles with reactions
   void SetNPartReac(Int_t np) {fNpReac=np;}
   void SetNPartCharge(Int_t nc) {fNpCharge=nc;}
   Int_t NPartReac() const {return fNpReac;}
   Int_t NPartCharge() const {return fNpCharge;}
   TDatabasePDG *DBPdg() const {return fDBPdg;}

   void SetEnergyGrid(Double_t emin, Double_t emax, Int_t nbins);
   Int_t    NEbins() const {return fNEbins;}
   Double_t Emin() const {return fEGrid[0];}
   Double_t Emax() const {return fEGrid[fNEbins-1];}
   Double_t EilDelta() const {return fEilDelta;}
   const Double_t* EGrid() const {return fEGrid;}

   static const char* EleSymb(Int_t z) {return fgEleSymbol[z-1];}
   static const char* EleName(Int_t z) {return fgEleName[z-1];}
   static Float_t WEle(Int_t z) {return fgWElem[z-1];}
   static Int_t NElem() {return fgNElem;}

   void Print(Option_t *option="") const;
   //approximated formula for nuclear mass computation; for handling fragments
   Double_t GetAprxNuclearMass(Int_t Z, Int_t A);

   void SetPDGToGVMap(std::map<Int_t,Int_t> &theMap);

   // only for e-,e+,gamma and proton 
   Int_t GetSpecGVIndex(Int_t indx){ return fSpecGVIndices[indx];}
   
private:
   TPartIndex(const TPartIndex&); // Not implemented
   TPartIndex& operator=(const TPartIndex&); // Not implemented
   
   static TPartIndex *fgPartIndex;

   const Int_t   fVersion=1000002;

   static const Int_t   fgNProc=FNPROC;    // Number of processes
   static const char   *fgPrName[FNPROC];  // Process name
   static const Short_t fgPCode[FNPROC];   // G4 process codes

   static const Int_t   fgNElem=NELEM;       // Number of Elements
   static const Char_t *fgEleSymbol[NELEM]; // Symbol of Element
   static const Char_t *fgEleName[NELEM];   // Name of Element
   static const Float_t fgWElem[NELEM];     // Weight of a mole in grams

   Int_t    fNPart;         // Total number of particles
   Int_t   *fPDG;           // [fNPart] PDG code of all part
   Int_t fNpReac;           // Number of particles with reactions
   Int_t fNpCharge;         // Number of particles with reactions

   Int_t     fNEbins;       // number of bins of common energy grid
   Double_t  fEilDelta;     // Inverse log delta of common energy grid
   Double_t *fEGrid;        // [fNEbins] Common energy grid

   TDatabasePDG *fDBPdg;    // Pointer to the augmented pdg database
   
   std::map<Int_t,Int_t> fPDGToGVMap;   // PDG->GV code map
   Int_t fSpecGVIndices[4]; // store GV codes of e-,e+,gamma and proton
   
   ClassDef(TPartIndex,1)  // Particle Index

};


#endif
