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
#include "TPFstate.h"
#include "Geant/Config.h"
#include "Geant/Error.h"

#ifndef GEANT_NVCC
#ifdef USE_ROOT
class TFile;
#endif
#endif

class TFinState;
class TPDecay;

#ifdef GEANT_NVCC
class TEFstate;
extern  int fEFNLdElemsHost;            //! number of loaded elements
extern GEANT_CUDA_DEVICE_CODE int fEFNLdElemsDev;            //! number of loaded elements
extern  TEFstate *fEFElementsHost[NELEM]; //! databases of elements
extern GEANT_CUDA_DEVICE_CODE TEFstate *fEFElementsDev[NELEM]; //! databases of elements
extern TPDecay  *fDecayHost;           //! decay table
extern GEANT_CUDA_DEVICE_CODE TPDecay  *fDecayDev;           //! decay table
#endif
class TEFstate {
public:
GEANT_CUDA_BOTH_CODE
  TEFstate();
  TEFstate(int z, int a, float dens);
  TEFstate &operator=(const TEFstate &other);
  ~TEFstate();
  static const char *ClassName() { return "TEFstate"; }
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[]);
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);

  bool AddPartFS(int kpart, int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                 const float en[], const char surv[], const int pid[], const float mom[]);
  GEANT_CUDA_BOTH_CODE
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

  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx) const;
  GEANT_CUDA_BOTH_CODE
  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom, double randn) const;
  GEANT_CUDA_BOTH_CODE
  bool GetReac(int pindex, int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
               const int *&pid, const float *&mom) const;

  static bool FloatDiff(double a, double b, double prec) { return fabs(a - b) > 0.5 * fabs(a + b) * prec; }
#ifndef GEANT_NVCC
  void Draw(const char *option);
#endif
  bool Resample();

  bool Prune();

GEANT_CUDA_BOTH_CODE
  int SizeOf() const;
  void Compact();
GEANT_CUDA_BOTH_CODE
  void RebuildClass();
GEANT_CUDA_BOTH_CODE
  static int SizeOfStore();
  static int MakeCompactBuffer(char* &b);
GEANT_CUDA_BOTH_CODE
  static void RebuildStore(char *b);
#ifdef MAGIC_DEBUG
GEANT_CUDA_BOTH_CODE
  int GetMagic() const {return fMagic;}
#endif
#ifndef GEANT_NVCC
  static int NLdElems() { return fNLdElems; }

  static TEFstate *Element(int i) {
    if (i < 0 || i >= fNLdElems)
      return 0;
    return fElements[i];
  }
#else
#ifdef GEANT_CUDA_DEVICE_BUILD
 GEANT_CUDA_DEVICE_CODE 
 static int NLdElems() { return fEFNLdElemsDev; }
 GEANT_CUDA_DEVICE_CODE 
 static TEFstate *Element(int i) {
    if (i < 0 || i >= fEFNLdElemsDev)
      return 0;
    return fEFElementsDev[i];
  }
#else
  static int NLdElems() { return fEFNLdElemsHost; }
  static TEFstate *Element(int i) {
    if (i < 0 || i >= fEFNLdElemsHost)
      return 0;
    return fEFElementsHost[i];
  }
#endif
#endif

#ifdef GEANT_NVCC
  GEANT_CUDA_BOTH_CODE
  static TEFstate *GetElement(int z, int a = 0);
#ifdef GEANT_DEVICE_BUILD
  GEANT_CUDA_DEVICE_CODE TEXsec **GetElements() { return fEFElementsDev; }
#else
  TEFstate **GetElements() { return fEFElementsHost; }
#endif
#else
  static TEFstate *GetElement(int z, int a = 0, TFile *f = 0);
  static TEFstate **GetElements() { return fElements; }
#endif
GEANT_CUDA_BOTH_CODE
  static TPDecay* GetDecayTable() {
#ifndef GEANT_NVCC
    return fDecay;
#else
#ifdef GEANT_CUDA_DEVICE_BUILD
    return fDecayDev;
#else
    return fDecayHost;
#endif
#endif
  }
GEANT_CUDA_BOTH_CODE
  static void SetDecayTable(TPDecay *decayTable) {
#ifndef GEANT_NVCC
    fDecay = decayTable;
#else
#ifdef GEANT_CUDA_DEVICE_BUILD
    fDecayDev = decayTable;
#else
    fDecayHost = decayTable;
#endif
#endif
  }

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
  TPFstate **fPFstateP; // [fNRpart] Final state table per particle

#ifndef GEANT_NVCC
  static int fNLdElems;              //! number of loaded elements
  static TEFstate *fElements[NELEM]; //! databases of elements
  static TPDecay  *fDecay;           //! decay table
#endif
#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
GEANT_CUDA_BOTH_CODE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) fEGrid) % sizeof(fEGrid[0]) != 0) { Geant::Error("TEFstate::CheckAlign","fEGrid misaligned\n");isaligned=false;}
  if(((unsigned long) &fAtcm3) % sizeof(fAtcm3) != 0) { Geant::Error("TEFstate::CheckAlign","fAtcm3 misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) { Geant::Error("TEFstate::CheckAlign","fEmin misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) { Geant::Error("TEFstate::CheckAlign","fEmax misaligned\n");isaligned=false;}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) { Geant::Error("TEFstate::CheckAlign","fEilDelta misaligned\n");isaligned=false;}
  if(((unsigned long) &fDens) % sizeof(fDens) != 0) { Geant::Error("TEFstate::CheckAlign","fDens misaligned\n");isaligned=false;}
  if(((unsigned long) &fEle) % sizeof(fEle) != 0) { Geant::Error("TEFstate::CheckAlign","fEle misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) { Geant::Error("TEFstate::CheckAlign","fNEbins misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEFstat) % sizeof(fNEFstat) != 0) { Geant::Error("TEFstate::CheckAlign","fNEFstat misaligned\n");isaligned=false;}
  if(((unsigned long) &fNRpart) % sizeof(fNRpart) != 0) { Geant::Error("TEFstate::CheckAlign","fNRpart misaligned\n");isaligned=false;}
  for(auto i=0; i< fNRpart; ++i)
    if(((unsigned long) fPFstateP[i]) % sizeof(double) != 0) { Geant::Error("TEFstate::CheckAlign","fPFstateP[%d misaligned\n",i);isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TEFstate::CheckAlign","fMagic misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TEFstate::CheckAlign","fStore misaligned\n");isaligned=false;}
  return isaligned;
}

#ifndef GEANT_NVCC
#ifdef USE_ROOT
  ClassDefNV(TEFstate, 4) // Element X-secs
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1];   // Pointer to the compact store part of the class
};

#endif
