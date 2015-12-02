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

#include "TPartIndex.h"
#ifdef USE_ROOT
#include "Rtypes.h"
#endif


class TFinState;

class TPFstate {
public:
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

  const char *Name() const { return TPartIndex::I()->PartName(fPDG); }
  bool SetPart(int pdg, int nfstat, int nreac, const int dict[]);
  bool SetPart(int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);
  bool SetFinState(int ibin, int reac, const int npart[], const float weight[], const float kerma[], const float en[],
                   const char surv[], const int pid[], const float mom[]);
  void Print(const char *opt = "") const;
  bool Prune() { return true; }
  bool SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx) const;
  bool SampleReac(int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom) const;
  bool SampleRestCaptFstate(int &npart, float &weight, float &kerma, float &enr, const int *&pid, const float *&mom,
                            double randn) const;

  bool GetReac(int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
               const float *&mom) const;
  int NEFstat() const { return fNEFstat; }
  void Dump() const {}
  bool Resample();

  int SizeOf() const;
  void Compact();
  void RebuildClass();
#ifdef MAGIC_DEBUG
  int GetMagic() const {return fMagic;}
#endif

bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) {std::cout << "TPFstate::fNEbins misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fNEFstat) % sizeof(fNEFstat) != 0) {std::cout << "TPFstate::fNEFstat misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fNFstat) % sizeof(fNFstat) != 0) {std::cout << "TPFstate::fNFstat misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fNReac) % sizeof(fNReac) != 0) {std::cout << "TPFstate::fNReac misaligned" << std::endl;isaligned=false;}
  for(auto i=0; i< fNFstat; ++i)
    if(((unsigned long) fFstatP[i]) % sizeof(double) != 0) {std::cout << "TPFstate::fFstatP[" << i << "] misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) fRestCaptFstat) % sizeof(double) != 0) {std::cout << "TPFstate::fRestCaptFstat misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fEGrid) % sizeof(fEGrid) != 0) {std::cout << "TPFstate::fEGrid misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) {std::cout << "TPFstate::fEmin misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) {std::cout << "TPFstate::fEmax misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) {std::cout << "TPFstate::fEilDelta misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fPDG) % sizeof(fPDG) != 0) {std::cout << "TPFstate::fPDG misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fRdict) % sizeof(int) != 0) {std::cout << "TPFstate::fRdict misaligned" << std::endl;isaligned=false;}
  if(((unsigned long) &fRmap) % sizeof(int) != 0) {std::cout << "TPFstate::fRmap misaligned" << std::endl;isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) {std::cout << "TPFstate::fMagic misaligned" << std::endl;isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) {std::cout << "TPFstate::fStore misaligned" << std::endl;isaligned=false;}
  return isaligned;
}

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }

private:
  TPFstate &operator=(const TPFstate &); // Not implemented

  static int fVerbose; // Controls verbosity level

  int fNEbins;               // number of energy bins
  int fNEFstat;              // number of states to sample per energy bin
  int fNFstat;               // tot size of fFstat
  int fNReac;                // number of reactions
  TFinState *fFstat;         // [fNFstat] table of final states
  TFinState **fFstatP;       //![fNFstat] table of pointers to final states
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
#ifdef USE_ROOT
  ClassDefNV(TPFstate, 2) // Particle Final States
#endif

private:
  alignas(sizeof(double)) char fStore[1]; // Pointer to compact memory
};

#endif
