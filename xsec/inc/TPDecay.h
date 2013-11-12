// @(#)root/base:$Id: $
// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TPDecay
#define ROOT_TPDecay


//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPDecay                                                              //
//                                                                      //
// Decay sampling for all particles                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Rtypes.h"

class TFinState;

class TPDecay {
public:
  
  TPDecay();
  TPDecay(Int_t nsample, Int_t npart, TFinState *decay);
  ~TPDecay();

  Int_t NSample() const {return fNSamp;}
  
  Bool_t SampleDecay(Int_t pindex, Int_t &npart, const Int_t *&pid, const Float_t *&mom) const;
  Bool_t GetDecay(Int_t pindex, Int_t ifs, Int_t &npart, const Int_t *&pid, const Float_t *&mom) const;
  
private:
  TPDecay(const TPDecay &); // Not implemented
  TPDecay& operator=(const TPDecay &); // Not implemented
  
  Int_t          fNSamp;   // Number of samples
  Int_t          fNPart;   // Number of particles
  TFinState     *fDecay;   //[fNPart] array of particle final states to be sampled
  
  ClassDefNV(TPDecay,1)  // Element X-secs
  
};


#endif
