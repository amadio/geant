// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TFinState
#define ROOT_TFinState


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TFinState                                                            //
//                                                                      //
// Final states for a reaction                                          //
//                                                                      //
// This class contains the final states for a given particle, energy    //
// and reaction                                                         //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TPartIndex.h"

class TFinState {
public:
  TFinState();
  TFinState(Int_t nfstates, const Int_t npart[], const Float_t weight[], const Float_t kerma[],
            const Float_t en[], const Char_t surv[], const Int_t pid[], const Float_t mom[]);
  ~TFinState();
  TFinState& operator=(const TFinState& right);
  
  Bool_t SetFinState(Int_t nfstates, const Int_t npart[], const Float_t weight[], const Float_t kerma[],
                     const Float_t en[], const Char_t surv[], const Int_t pid[], const Float_t mom[]);
  Bool_t SetFinState(const TFinState &right);
  void NormFinSateWeights();
  Int_t GetNsecs() const {return fNsecs;}

  Bool_t Prune() {return kTRUE;}
  Bool_t SampleReac(Int_t& npart, Float_t& weight, Float_t& kerma, Float_t &en, const Int_t *&pid,
                    const Float_t *&mom) const;
  Bool_t SampleReac(Int_t& npart, Float_t& weight, Float_t& kerma, Float_t &en, const Int_t *&pid,
                    const Float_t *&mom, Double_t randn) const;

  Bool_t GetReac(Int_t finstat, Int_t& npart, Float_t& weight, Float_t& kerma, Float_t &en,
                 const Int_t *&pid, const Float_t *&mom) const;
  void Dump() const {}
  void Print(Option_t */*opt*/="") const
  { printf("fNFstates %d, fNsecs %d, fNMom %d, fPID %p, fSurv %p, fNpart %p, fWeight %p, fKerma %p, fMom %p\n",
           fNFstates,fNsecs,fNMom,(void *)fPID,(void *)fSurv,
           (void *)fNpart,(void *)fWeight,(void *)fKerma, (void *)fMom);}
  
  static void SetVerbose(Int_t verbose) {fVerbose=verbose;}
  static Int_t GetVerbose() {return fVerbose;}
  
private:
  TFinState(const TFinState&);    // Not implemented
  
  static Int_t    fVerbose;       // Controls verbosity level
  
  Int_t           fNFstates;      // Number of final states
  Int_t           fNsecs;         // Total number of secondaries
  Int_t           fNMom;          // 3*fNsecs, just because ROOT cannot use formulas in dimensions
  Int_t          *fNpart;         // [fNFstates] number of particles in each final state
  Float_t        *fWeight;        // [fNFstates] Weight of the final states
  Float_t        *fKerma;         // [fNFstates] Released energy
  Float_t        *fEn;            // [fNFstates] Energy of final states in GeV
  Char_t         *fSurv;          // [fNFstates] whether the orignal particle has survived or not
  Int_t          *fPID;           // [fNsecs] GeantV particle code
  Float_t        *fMom;           // [fNMom] Particle momentum (GeV)
  
  ClassDefNV(TFinState,1)  //Particle Final States
};

#endif
