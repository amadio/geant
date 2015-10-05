// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TPDecay_H
#define TPDecay_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPDecay                                                              //
//                                                                      //
// Decay sampling for all particles                                     //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef USE_ROOT
#include "RTypes.h"
#endif

class TFinState;

class TPDecay {
public:
  TPDecay();
  TPDecay(int nsample, int npart, TFinState *decay);
  ~TPDecay();

  int NSample() const { return fNSamp; }

  bool SampleDecay(int pindex, int &npart, const int *&pid, const float *&mom) const;
  bool GetDecay(int pindex, int ifs, int &npart, const int *&pid, const float *&mom) const;
  bool HasDecay(int pindex) const;
  void SetCTauPerMass(double *ctaupermass, int np);
  double GetCTauPerMass(int pindex) const {
    // should check if it is initialized but we know what we are doing!
    return fCTauPerMass[pindex];
  }

private:
  TPDecay(const TPDecay &);            // Not implemented
  TPDecay &operator=(const TPDecay &); // Not implemented

  int fNSamp;           // Number of samples
  int fNPart;           // Number of particles
  TFinState *fDecay;    // [fNPart] array of particle final states to be sampled
  double *fCTauPerMass; // [fNPart] precomputed c*tau/mass values [cm/GeV]

#ifdef USE_ROOT
  ClassDefNV(TPDecay, 1) // Element X-secs
#endif
};

#endif
