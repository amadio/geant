// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TPFstate_H
#define TPFstate_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPFstate                                                             //
//                                                                      //
// Final states for the reactions of a particle                         //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "Geant/Config.h"
#include "TPartIndex.h"
#include "Geant/Error.h"

#ifndef GEANT_NVCC
#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#endif

class TFinState;

class TPFstate {
public:
  GEANT_CUDA_BOTH_CODE
  TPFstate();
  TPFstate(int pdg, int nfstat, int nreac, const int dict[]);
  TPFstate(const TPFstate &other);
  ~TPFstate();

  void SetRestCaptFstate(const TFinState &finstate);
  bool HasRestCaptFstat() {
    if (!fRestCaptFstat)
      return false;
    return true;
  }
  GEANT_CUDA_BOTH_CODE
  const char *Name() const { return TPartIndex::I()->PartName(fPDG); }

  bool SetPart(int pdg, int nfstat, int nreac, const int dict[]);
  bool SetPart(int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);
  bool SetFinState(int ibin, int reac, const int npart[], const float weight[], const float kerma[], const float en[],
                   const char surv[], const int pid[], const float mom[]);
  void Print(const char *opt = "") const;
  bool Prune() { return true; }
  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx) const;
  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom) const;
  bool SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid, const float *&mom,
                            double randn) const;

  GEANT_CUDA_BOTH_CODE
  bool GetReac(int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
               const float *&mom) const;
  int NEFstat() const { return fNEFstat; }
  void Dump() const {}
  bool Resample();

  GEANT_CUDA_BOTH_CODE
  int SizeOf() const;
  void Compact();
  GEANT_CUDA_BOTH_CODE
  void RebuildClass();
#ifdef MAGIC_DEBUG
  GEANT_CUDA_BOTH_CODE
  int GetMagic() const { return fMagic;}
#endif

  GEANT_CUDA_BOTH_CODE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) { Geant::Error("TPFstate::CheckAlign","fNEbins misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEFstat) % sizeof(fNEFstat) != 0) { Geant::Error("TPFstate::CheckAlign","fNEFstat misaligned\n");isaligned=false;}
  if(((unsigned long) &fNFstat) % sizeof(fNFstat) != 0) { Geant::Error("TPFstate::CheckAlign","fNFstat misaligned\n");isaligned=false;}
  if(((unsigned long) &fNReac) % sizeof(fNReac) != 0) { Geant::Error("TPFstate::CheckAlign","fNReac misaligned\n");isaligned=false;}
  for(auto i=0; i< fNFstat; ++i)
    if(((unsigned long) fFstatP[i]) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","fFstatP[%d] misaligned\n",i);isaligned=false;}
  if(((unsigned long) fRestCaptFstat) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","fRestCaptFstat misaligned\n");isaligned=false;}
  if(((unsigned long) &fEGrid) % sizeof(fEGrid) != 0) { Geant::Error("TPFstate::CheckAlign","fEGrid misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) { Geant::Error("TPFstate::CheckAlign","fEmin misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) { Geant::Error("TPFstate::CheckAlign","fEmax misaligned\n");isaligned=false;}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) { Geant::Error("TPFstate::CheckAlign","fEilDelta misaligned\n");isaligned=false;}
  if(((unsigned long) &fPDG) % sizeof(fPDG) != 0) { Geant::Error("TPFstate::CheckAlign","fPDG misaligned\n");isaligned=false;}
  if(((unsigned long) &fRdict) % sizeof(int) != 0) { Geant::Error("TPFstate::CheckAlign","fRdict misaligned\n");isaligned=false;}
  if(((unsigned long) &fRmap) % sizeof(int) != 0) { Geant::Error("TPFstate::CheckAlign","fRmap misaligned\n");isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TPFstate::CheckAlign","fMagic misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","fStore misaligned\n");isaligned=false;}
  return isaligned;
}

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }

private:
  TPFstate &operator=(const TPFstate &) = delete;

  static int fVerbose; // Controls verbosity level

  int fNEbins;               // number of energy bins
  int fNEFstat;              // number of states to sample per energy bin
  int fNFstat;               // tot size of fFstat
  int fNReac;                // number of reactions
  TFinState *fFstat;         // [fNFstat] table of final states
  TFinState **fFstatP;       // [fNFstat] table of pointers to final states
  TFinState *fRestCaptFstat; // RestCapture final states
  const double *fEGrid;      //![fNEbins] energy grid
  double fEmin;              // Min energy of the energy grid
  double fEmax;              // Max energy of the energy grid
  double fEilDelta;          // logarithmic energy delta
  int fPDG;                  // particle pdg code

  int fRdict[FNPROC]; // reaction dictionary from reaction number to position
  // in the X-sec array
  int fRmap[FNPROC]; // reaction map, from reaction position in the X-sec
// array to the raction number

#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
#ifndef GEANT_NVCC
#ifdef USE_ROOT
  ClassDefNV(TPFstate, 3) // Particle Final States
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1]; // Pointer to compact memory
};

#endif
