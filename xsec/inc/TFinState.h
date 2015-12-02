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
#ifdef USE_ROOT
#include "Rtypes.h"
#endif

class TFinState {
public:
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
  bool SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid, const float *&mom) const;
  bool SampleReac(int &npart, float &weight, float &kerma, float &en, const int *&pid, const float *&mom,
                  double randn) const;

  bool GetReac(int finstat, int &npart, float &weight, float &kerma, float &en, const int *&pid,
               const float *&mom) const;
  void Dump() const {}
  void Print(const char * /*opt*/ = "") const {
    printf("fNFstates %d, fNsecs %d, fNMom %d, fPID %p, fSurv %p, fNpart %p, fWeight %p, fKerma %p, fMom %p\n",
           fNFstates, fNsecs, fNMom, (void *)fPID, (void *)fSurv, (void *)fNpart, (void *)fWeight, (void *)fKerma,
           (void *)fMom);
  }

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }

  int SizeOf() const;
  void Compact();
  void RebuildClass();
#ifdef MAGIC_DEBUG
  int GetMagic() const {return fMagic;}
#endif

bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fNFstates) % sizeof(fNFstates) != 0) {std::cout << "TPFstate::fNFstates misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fNsecs) % sizeof(fNsecs) != 0) {std::cout << "TPFstate::fNsecs misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fNMom) % sizeof(fNMom) != 0) {std::cout << "TPFstate::fNMom misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fWeight) % sizeof(fWeight[0]) != 0) {std::cout << "TPFstate::fWeight misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fKerma) % sizeof(fKerma[0]) != 0) {std::cout << "TPFstate::fKerma misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fEn) % sizeof(fEn[0]) != 0) {std::cout << "TPFstate::fEn misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fMom) % sizeof(fMom[0]) != 0) {std::cout << "TPFstate::fMom misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fPID) % sizeof(fPID[0]) != 0) {std::cout << "TPFstate::fPID misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fNpart) % sizeof(fNpart[0]) != 0) {std::cout << "TPFstate::fNpart misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fSurv) % sizeof(fSurv[0]) != 0) {std::cout << "TPFstate::fSurv misaligned" << std::endl;isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) {std::cout << "TPFstate::fMagic misaligned" << std::endl;isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) {std::cout << "TPFstate::fStore misaligned" << std::endl;isaligned=false;}
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
#ifdef USE_ROOT
  ClassDefNV(TFinState, 2) // Particle Final States
#endif

private:
  alignas(sizeof(double)) char fStore[1];        //! Pointer to the compact data of the class
};

#endif
