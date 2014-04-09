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


#include "TDatabasePDG.h"
#include "TPartIndex.h"

class TFinState;

class TPFstate {
public:
  TPFstate();
  TPFstate(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]);
  ~TPFstate();
  
  void SetRestCaptFstate(const TFinState &finstate);

  const char* Name() const {return TDatabasePDG::Instance()->GetParticle(fPDG)->GetName();}
  Bool_t SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[]);
  Bool_t SetPart(Int_t pdg, Int_t nfstat, Int_t nreac, const Int_t dict[], TFinState vecfs[]);
  Bool_t SetFinState(Int_t ibin, Int_t reac, const Int_t npart[], const Float_t weight[], const Float_t kerma[],
                     const Float_t en[], const Char_t surv[], const Int_t pid[], const Float_t mom[]);
  void Print(Option_t *opt="") const;
  Bool_t Prune() {return kTRUE;}
  Bool_t SampleReac(Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const;
  Bool_t SampleReac(Int_t preac, Float_t en, Int_t& npart, Float_t& weight,
                    Float_t& kerma, Float_t &enr, const Int_t *&pid,
                    const Float_t *&mom, Double_t randn1, Double_t randn2) const;
  Bool_t SampleRestCaptFstate(Int_t& npart, Float_t& weight, Float_t& kerma, 
		    Float_t &enr, const Int_t *&pid, const Float_t *&mom) const;
  Bool_t SampleRestCaptFstate(Int_t& npart, Float_t& weight, Float_t& kerma,
                    Float_t &enr, const Int_t *&pid, const Float_t *&mom,
                    Double_t randn) const;

  Bool_t GetReac(Int_t preac, Float_t en, Int_t ifs, Int_t& npart, Float_t& weight,
                 Float_t& kerma, Float_t &enr, const Int_t *&pid, const Float_t *&mom) const;
  Int_t NEFstat() const {return fNEFstat;}
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
  TFinState      *fRestCaptFstat; // RestCapture final states

  Int_t           fRdict[FNPROC]; // reaction dictionary from reaction number to position
  // in the X-sec array
  Int_t           fRmap[FNPROC];  // reaction map, from reaction position in the X-sec
  // array to the raction number
  
  ClassDefNV(TPFstate,1)  //Particle Final States
};

#endif

