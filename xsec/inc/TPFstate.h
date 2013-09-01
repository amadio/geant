// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPFstate
#define ROOT_TPFstate


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPFstate                                                             //
//                                                                      //
// Final states for the reactions of a particle                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include <TObject.h>
#include <TDatabasePDG.h>
#include <TPartIndex.h>
class TFinState;

class TPFstate: public TObject {
public:
  TPFstate();
  TPFstate(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]);
  ~TPFstate();
  
  const char* Name() const {return TDatabasePDG::Instance()->GetParticle(fPDG)->GetName();}
  Bool_t SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]);
  Bool_t SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[], TFinState vecfs[]);
  Bool_t SetFinState(Double_t en, Int_t reac, const Float_t weight[], const Float_t kerma[],
                     const Int_t npart[], const Float_t mom[], const Int_t pid[], const Char_t surv[]);
  void Print(Option_t *opt="") const;
  Bool_t Prune() {return kTRUE;}
  Bool_t SampleReac(Double_t en, Int_t preac, Float_t& kerma, Int_t& npart, Int_t *pid, Float_t mom[]) const;
  void Dump() const {}
  Bool_t Resample();
  
  static void SetVerbose(Int_t verbose) {fVerbose=verbose;}
  static Int_t GetVerbose() {return fVerbose;}
  
private:
  TPFstate(const TPFstate&);      // Not implemented
  TPFstate& operator=(const TPFstate&);      // Not implemented
  
  static Int_t    fVerbose;       // Controls verbosity level
  
  Int_t           fPDG;           // particle pdg code
  Int_t           fNEbins;        // number of energy bins
  Int_t           fNReac;         // number of reactions
  Int_t           fNEFstat;       // number of states to sample per energy bin
  Int_t           fNFstat;        // tot size of fFstat
  Double_t        fEmin;          // Min energy of the energy grid
  Double_t        fEmax;          // Max energy of the energy grid
  Double_t        fEilDelta;      // logarithmic energy delta
  const Double_t *fEGrid;         //![fNEbins] energy grid
  TFinState      *fFstat;         // [fNFstat] table of final states
  Int_t           fRdict[FNPROC]; // reaction dictionary from reaction number to position
  // in the X-sec array
  Int_t           fRmap[FNPROC];  // reaction map, from reaction position in the X-sec
  // array to the raction number
  
  ClassDef(TPFstate,1)  //Particle Final States
};

#endif

