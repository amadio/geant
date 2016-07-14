// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TFinState_H
#define TFinState_H

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
#include "Geant/Error.h"

class TFinState {
public: 
 GEANT_CUDA_BOTH_CODE
  TFinState();
  TFinState(int nfstates, const int npart[], const float weight[], const float kerma[], const float en[],
            const char surv[], const int pid[], const float mom[]);
  TFinState(const TFinState & other);
  ~TFinState();
  TFinState &operator=(const TFinState &right);

  bool SetFinState(int nfstates, const int npart[], const float weight[], const float kerma[], const float en[],
                   const char surv[], const int pid[], const float mom[]);
  bool SetFinState(const TFinState &right);
  void NormFinSateWeights();
  int GetNsecs() const { return fNsecs; }

  bool Prune() { return true; }
  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid, const float *&mom) const;
  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid, const float *&mom,
                  double randn) const;

  GEANT_CUDA_BOTH_CODE
  bool GetReac(int finstat, int &npart, float &weight, float &kerma, float &en, const int *&pid,
               const float *&mom) const;
  void Dump() const {}
  void Print(const char * /*opt*/ = "") const {
    Geant::Printf("fNFstates %d, fNsecs %d, fNMom %d, fPID %p, fSurv %p, fNpart %p, fWeight %p, fKerma %p, fMom %p\n",
                  fNFstates, fNsecs, fNMom, (void *)fPID, (void *)fSurv, (void *)fNpart, (void *)fWeight, (void *)fKerma,
                  (void *)fMom);
  }

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }

  GEANT_CUDA_BOTH_CODE
  int SizeOf() const;
  void Compact();
  GEANT_CUDA_BOTH_CODE
  void RebuildClass();
#ifdef MAGIC_DEBUG
  GEANT_CUDA_BOTH_CODE
  int GetMagic() const {return fMagic;}
#endif

  GEANT_CUDA_BOTH_CODE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fNFstates) % sizeof(fNFstates) != 0) { Geant::Error("TPFstate::CheckAlign","fNFstates misaligned\n");isaligned=false;}
  if(((unsigned long) &fNsecs) % sizeof(fNsecs) != 0) { Geant::Error("TPFstate::CheckAlign","fNsecs misaligned\n");isaligned=false;}
  if(((unsigned long) &fNMom) % sizeof(fNMom) != 0) { Geant::Error("TPFstate::CheckAlign","fNMom misaligned\n");isaligned=false;}
  if(((unsigned long) fWeight) % sizeof(fWeight[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fWeight misaligned\n");isaligned=false;}
  if(((unsigned long) fKerma) % sizeof(fKerma[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fKerma misaligned\n");isaligned=false;}
  if(((unsigned long) fEn) % sizeof(fEn[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fEn misaligned\n");isaligned=false;}
  if(((unsigned long) fMom) % sizeof(fMom[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fMom misaligned\n");isaligned=false;}
  if(((unsigned long) fPID) % sizeof(fPID[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fPID misaligned\n");isaligned=false;}
  if(((unsigned long) fNpart) % sizeof(fNpart[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fNpart misaligned\n");isaligned=false;}
  if(((unsigned long) fSurv) % sizeof(fSurv[0]) != 0) { Geant::Error("TPFstate::CheckAlign","fSurv misaligned\n");isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TPFstate::CheckAlign","fMagic misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TPFstate::CheckAlign","fStore misaligned\n");isaligned=false;}
  return isaligned;
}

private:
  static int fVerbose; // Controls verbosity level

  int fNFstates;  // Number of final states
  int fNsecs;     // Total number of secondaries
  int fNMom;      // 3*fNMom, just because ROOT cannot use formulas in dimensions
  float *fWeight; // [fNFstates] Weight of the final states
  float *fKerma;  // [fNFstates] Released energy
  float *fEn;     // [fNFstates] Energy of final states in GeV
  float *fMom;    // [fNMom] Particle momentum (GeV)
  int *fPID;      // [fNsecs] GeantV particle code
  int *fNpart;    // [fNFstates] number of particles in each final state
  char *fSurv;    // [fNFstates] whether the orignal particle has survived or not

#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif

private:
  alignas(sizeof(double)) char fStore[1];        //! Pointer to the compact data of the class
};

#endif
