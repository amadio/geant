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
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#endif
#include "Geant/Config.h"
#include "Geant/Error.h"

#include <iostream>

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

  VECCORE_ATT_HOST_DEVICE
  int SizeOf() const;
  void Compact();
  VECCORE_ATT_HOST_DEVICE
  void RebuildClass();
  size_t MakeCompactBuffer(char* &b);
#ifdef MAGIC_DEBUG
  VECCORE_ATT_HOST_DEVICE
  int GetMagic() const {return fMagic;}
#endif

  VECCORE_ATT_HOST_DEVICE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fNPart) % sizeof(fNPart) != 0) { Geant::Error("TPFstate::CheckAlign","%s","fNPart misaligned\n");isaligned=false;}
  if(((unsigned long) fCTauPerMass) % sizeof(fCTauPerMass) != 0) { Geant::Error("TPFstate::CheckAlign","%s","fCTauPerMass misaligned\n");isaligned=false;}
  if(((unsigned long) fDecayP) % sizeof(fDecayP) != 0) { Geant::Error("TPFstate::CheckAlign","%s","fDecayP misaligned\n");isaligned=false;}
  for(auto i=0; i< fNPart; ++i)
    if(((unsigned long) fDecayP[i]) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","fDecayP[%d] misaligned\n",i);isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TPFstate::CheckAlign","%s","fMagic misaligned");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","%s","fStore misaligned");isaligned=false;}
  return isaligned;
}

private:
#ifndef VECCORE_CUDA
  TPDecay &operator=(const TPDecay &); // Not implemented
#endif

  int fNPart;           // Number of particles
  double *fCTauPerMass; // [fNPart] precomputed c*tau/mass values [cm/GeV]
  TFinState *fDecay;    // [fNPart] array of particle final states to be sampled
  TFinState **fDecayP;  // [fNPart] table of pointers to final states

#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
#ifndef VECCORE_CUDA
#ifdef USE_ROOT
  ClassDefNV(TPDecay, 5) // Element X-secs
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1];    // Pointer to the compact part of the store
};

#endif
