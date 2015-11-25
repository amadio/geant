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
#include "Rtypes.h"
#endif

class TFinState;

class TPDecay {
public:
  TPDecay();
  TPDecay(int nsample, int npart, TFinState *decay);
  TPDecay(const TPDecay &other);
  ~TPDecay();

  bool SampleDecay(int pindex, int &npart, const int *&pid, const float *&mom) const;
  bool GetDecay(int pindex, int ifs, int &npart, const int *&pid, const float *&mom) const;
  bool HasDecay(int pindex) const;
  void SetCTauPerMass(double *ctaupermass, int np);
  double GetCTauPerMass(int pindex) const {
    // should check if it is initialized but we know what we are doing!
    return fCTauPerMass[pindex];
  }

  int SizeOf() const;
  void Compact();
  void RebuildClass();
  size_t MakeCompactBuffer(char* &b);
#ifdef MAGIC_DEBUG
  int GetMagic() const {return fMagic;}
#endif

private:
  TPDecay &operator=(const TPDecay &); // Not implemented

  int fNPart;           // Number of particles
  double *fCTauPerMass; // [fNPart] precomputed c*tau/mass values [cm/GeV]
  TFinState *fDecay;    // [fNPart] array of particle final states to be sampled
  TFinState **fDecayP;  //![fNpart] table of pointers to final states

#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
#ifdef USE_ROOT
  ClassDefNV(TPDecay, 2) // Element X-secs
#endif
  double fStore[1];    // Pointer to the compact part of the store
};

#endif
