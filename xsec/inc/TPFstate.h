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

class TFinState;

class TPFstate {
public:
  TPFstate();
  TPFstate(int pdg, int nfstat, int nreac, const int dict[]);
  ~TPFstate();

  void SetRestCaptFstate(const TFinState &finstate);
  bool HasRestCaptFstat() {
    if (!fRestCaptFstat)
      return kFALSE;
    return kTRUE;
  }

#ifdef USE_VECGEOM_NAVIGATOR
  const char *Name() const { return Particle::GetParticle(fPDG).Name(); }
#else
  const char *Name() const { return TDatabasePDG::Instance()->GetParticle(fPDG)->GetName(); }
#endif
  bool SetPart(int pdg, int nfstat, int nreac, const int dict[]);
  bool SetPart(int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);
  bool SetFinState(int ibin, int reac, const int npart[], const Float_t weight[], const Float_t kerma[],
                   const Float_t en[], const Char_t surv[], const int pid[], const Float_t mom[]);
  void Print(Option_t *opt = "") const;
  bool Prune() { return kTRUE; }
  bool SampleReac(int preac, Float_t en, int &npart, Float_t &weight, Float_t &kerma, Float_t &enr, const int *&pid,
                  const Float_t *&mom, int &ebinindx) const;
  bool SampleReac(int preac, Float_t en, int &npart, Float_t &weight, Float_t &kerma, Float_t &enr, const int *&pid,
                  const Float_t *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int &npart, Float_t &weight, Float_t &kerma, Float_t &enr, const int *&pid,
                            const Float_t *&mom) const;
  bool SampleRestCaptFstate(int &npart, Float_t &weight, Float_t &kerma, Float_t &enr, const int *&pid,
                            const Float_t *&mom, double randn) const;

  bool GetReac(int preac, Float_t en, int ifs, int &npart, Float_t &weight, Float_t &kerma, Float_t &enr,
               const int *&pid, const Float_t *&mom) const;
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

  ClassDefNV(TPFstate, 1) // Particle Final States
};

#endif
