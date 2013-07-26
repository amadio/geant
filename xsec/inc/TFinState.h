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


#include <TObject.h>
#include <TDatabasePDG.h>
#include <TPartIndex.h>

class TFinState: public TObject {
public:
   TFinState();
   TFinState(Int_t nfstates, const Float_t weight[],
	     const Float_t kerma[], const Int_t npart[],
	     const Float_t (*mom)[3], const Int_t pid[],
             const Char_t surv[]);
   ~TFinState();
   TFinState& operator=(const TFinState& right);

   Bool_t SetFinState(Int_t nfstates, const Float_t weight[],
	     const Float_t kerma[], const Int_t npart[],
	     const Float_t (*mom)[3], const Int_t pid[],
                      const Char_t surv[]);
   void Print(Option_t */*opt*/="") const {}
   Bool_t Prune() {return kTRUE;}
   Bool_t SampleReac(Int_t& npart, Int_t* pid, Float_t (*mom)[3]) const;
   void Dump() const {}

   static void SetVerbose(Int_t verbose) {fVerbose=verbose;}
   static Int_t GetVerbose() {return fVerbose;}

private:
   TFinState(const TFinState&);    // Not implemented

   static Int_t    fVerbose;       // Controls verbosity level

   Int_t           fNFstates;      // Number of final states
   Int_t           fNsecs;         // Total number of secondaries
   Int_t          *fPID;           // [fNsecs] G5 particle code
   Char_t         *fSurv;          // [fNFstates] whether the orignal particle has survived or not
   Int_t          *fNpart;         // [fNFstates] number of particles in each final state
   Float_t        *fWeight;        // [fNFstates] Weight of the final states
   Float_t        *fKerma;         // [fNFstates] Released energy
   Float_t      (*fMom)[3];        // [fNsecs] Particle momentum (GeV)

   ClassDef(TFinState,1)  //Particle Final States
};

#endif
