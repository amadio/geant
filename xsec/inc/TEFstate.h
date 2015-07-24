// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TEFstate_H
#define TEFstate_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TEXSec                                                               //
//                                                                      //
// X-section for GV per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPartIndex.h"

class TFile;

class TPFstate;
class TFinState;

class TEFstate {
public:
  TEFstate();
  TEFstate(int z, int a, float dens);
  ~TEFstate();
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[]);
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);

  bool AddPartFS(int kpart, int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                 const float en[], const char surv[], const int pid[], const float mom[]);

  int Ele() const { return fEle; }
  double Dens() const { return fDens; }
  double Emin() const { return fEmin; }
  double Emax() const { return fEmax; }
  int NEbins() const { return fNEbins; }
  double EilDelta() const { return fEilDelta; }
  int NEFstat() const { return fNEFstat; }

  int NRpart() const { return fNRpart; }

  void SetRestCaptFstate(int kpart, const TFinState &fstate);
  bool HasRestCapture(int partindex);

  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx) const;
  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom, double randn) const;
  bool GetReac(int pindex, int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
               const int *&pid, const float *&mom) const;

  static bool FloatDiff(double a, double b, double prec) { return fabs(a - b) > 0.5 * fabs(a + b) * prec; }

  void Draw(const char *option);
  bool Resample();

  bool Prune();

  static int NLdElems() { return fNLdElems; }

  static TEFstate *Element(int i) {
    if (i < 0 || i >= fNLdElems)
      return 0;
    return fElements[i];
  }

  static TEFstate *GetElement(int z, int a = 0, TFile *f = 0);
  static TEFstate **GetElements() { return fElements; }

private:
  TEFstate(const TEFstate &);            // Not implemented
  TEFstate &operator=(const TEFstate &); // Not implemented

  int fEle;             // Element code Z*10000+A*10+metastable level
  float fDens;          // Density in g/cm3
  double fAtcm3;        // Atoms per cubic cm unit density
  double fEmin;         // Minimum of the energy Grid
  double fEmax;         // Maximum of the energy Grid
  int fNEbins;          // Number of log steps in energy
  double fEilDelta;     // Inverse log energy step
  const double *fEGrid; //! Common energy grid
  int fNEFstat;         // Number of sampled final states
  int fNRpart;          // Number of particles with reaction
  TPFstate *fPFstate;   // [fNRpart] Final state table per particle

  static int fNLdElems;              //! number of loaded elements
  static TEFstate *fElements[NELEM]; //! databases of elements

  ClassDefNV(TEFstate, 1) // Element X-secs
};

#endif
