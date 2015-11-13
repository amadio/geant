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

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }

private:
  TPFstate(const TPFstate &);            // Not implemented
  TPFstate &operator=(const TPFstate &); // Not implemented

  static int fVerbose; // Controls verbosity level

  int fPDG;                  // particle pdg code
  int fNEbins;               // number of energy bins
  int fNReac;                // number of reactions
  int fNEFstat;              // number of states to sample per energy bin
  int fNFstat;               // tot size of fFstat
  double fEmin;              // Min energy of the energy grid
  double fEmax;              // Max energy of the energy grid
  double fEilDelta;          // logarithmic energy delta
  const double *fEGrid;      //![fNEbins] energy grid
  TFinState *fFstat;         // [fNFstat] table of final states
  TFinState *fRestCaptFstat; // RestCapture final states

  int fRdict[FNPROC]; // reaction dictionary from reaction number to position
  // in the X-sec array
  int fRmap[FNPROC]; // reaction map, from reaction position in the X-sec
// array to the raction number

#ifdef USE_ROOT
  ClassDefNV(TPFstate, 1) // Particle Final States
#endif
};

#endif
