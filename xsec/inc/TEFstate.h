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
#include "TPFstate.h"

class TFile;

class TFinState;
class TPDecay;

class TEFstate {
public:
  TEFstate();
  TEFstate(int z, int a, float dens);
  TEFstate &operator=(const TEFstate &other);
  ~TEFstate();
  static const char *ClassName() { return "TEFstate"; }
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

  int SizeOf() const;
  void Compact();
  void RebuildClass();
  static int SizeOfStore();
  static int MakeCompactBuffer(char* &b);
  static void RebuildStore(char *b);
#ifdef MAGIC_DEBUG
  int GetMagic() const {return fMagic;}
#endif

  static int NLdElems() { return fNLdElems; }

  static TEFstate *Element(int i) {
    if (i < 0 || i >= fNLdElems)
      return 0;
    return fElements[i];
  }

  static TEFstate *GetElement(int z, int a = 0, TFile *f = 0);
  static TEFstate **GetElements() { return fElements; }
  static TPDecay* GetDecayTable() {return fDecay;}
  static void SetDecayTable(TPDecay *decayTable) {fDecay = decayTable;}

private:
  TEFstate(const TEFstate &);            // Not implemented

  const double *fEGrid; //! Common energy grid
  double fAtcm3;        // Atoms per cubic cm unit density
  double fEmin;         // Minimum of the energy Grid
  double fEmax;         // Maximum of the energy Grid
  double fEilDelta;     // Inverse log energy step
  float fDens;          // Density in g/cm3
  int fEle;             // Element code Z*10000+A*10+metastable level
  int fNEbins;          // Number of log steps in energy
  int fNEFstat;         // Number of sampled final states
  int fNRpart;          // Number of particles with reaction
  TPFstate *fPFstate;   // [fNRpart] Final state table per particle
  TPFstate **fPFstateP; //! [fNRpart] Final state table per particle

  static int fNLdElems;              //! number of loaded elements
  static TEFstate *fElements[NELEM]; //! databases of elements
  static TPDecay  *fDecay;           //! decay table
#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif

void CheckAlign() {
  if(((unsigned long) fEGrid) % sizeof(fEGrid[0]) != 0) {std::cout << "TEFstate::fEGrid misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fAtcm3) % sizeof(fAtcm3) != 0) {std::cout << "TEFstate::fAtcm3 misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) {std::cout << "TEFstate::fEmin misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) {std::cout << "TEFstate::fEmax misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) {std::cout << "TEFstate::fEilDelta misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fDens) % sizeof(fDens) != 0) {std::cout << "TEFstate::fDens misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fEle) % sizeof(fEle) != 0) {std::cout << "TEFstate::fEle misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) {std::cout << "TEFstate::fNEbins misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fNEFstat) % sizeof(fNEFstat) != 0) {std::cout << "TEFstate::fNEFstat misaligned" << std::endl;exit(1);}
  if(((unsigned long) &fNRpart) % sizeof(fNRpart) != 0) {std::cout << "TEFstate::fNRpart misaligned" << std::endl;exit(1);}
  for(auto i=0; i< fNRpart; ++i)
    if(((unsigned long) fPFstateP[i]) % sizeof(double) != 0) {std::cout << "TEFstate::fPFstateP[" << i << "] misaligned" << std::endl;exit(1);}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) {std::cout << "TEFstate::fMagic misaligned" << std::endl;exit(1);}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) {std::cout << "TEFstate::fStore misaligned" << std::endl;exit(1);}
}

#ifdef USE_ROOT
  ClassDefNV(TEFstate, 2) // Element X-secs
#endif

private:
  alignas(sizeof(double)) char fStore[1];   // Pointer to the compact store part of the class
};

#endif
