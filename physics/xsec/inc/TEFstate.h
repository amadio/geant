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

#include "Geant/Config.h"
#include "Geant/Error.h"
#include "TPFstate.h"

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
class TFile;
#endif
#endif

class TFinState;
class TPDecay;

#ifdef VECCORE_CUDA
class TEFstate;
extern  int fEFNLdElemsHost;            //! number of loaded elements
extern VECCORE_ATT_DEVICE int fEFNLdElemsDev;            //! number of loaded elements
extern  TEFstate *fEFElementsHost[NELEM]; //! databases of elements
extern VECCORE_ATT_DEVICE TEFstate *fEFElementsDev[NELEM]; //! databases of elements
extern TPDecay  *fDecayHost;           //! decay table
extern VECCORE_ATT_DEVICE TPDecay  *fDecayDev;           //! decay table
#endif

class TEFstate {
public:
  VECCORE_ATT_HOST_DEVICE
  TEFstate();
  TEFstate(int z, int a, float dens);
  TEFstate &operator=(const TEFstate &other);
  ~TEFstate();
  static const char *ClassName() { return "TEFstate"; }
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[]);
  bool AddPart(int kpart, int pdg, int nfstat, int nreac, const int dict[], TFinState vecfs[]);

  bool AddPartFS(int kpart, int ibin, int reac, const int npart[], const float weight[], const float kerma[],
                 const float en[], const char surv[], const int pid[], const float mom[]);
  VECCORE_ATT_HOST_DEVICE
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

  VECCORE_ATT_HOST_DEVICE
  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx) const;
  VECCORE_ATT_HOST_DEVICE
  bool SampleReac(int pindex, int preac, float en, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                  const float *&mom, int &ebinindx, double randn1, double randn2) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom) const;
  bool SampleRestCaptFstate(int kpart, int &npart, float &weight, float &kerma, float &enr, const int *&pid,
                            const float *&mom, double randn) const;
  VECCORE_ATT_HOST_DEVICE
  bool GetReac(int pindex, int preac, float en, int ifs, int &npart, float &weight, float &kerma, float &enr,
               const int *&pid, const float *&mom) const;

  static bool FloatDiff(double a, double b, double prec) { return fabs(a - b) > 0.5 * fabs(a + b) * prec; }
#ifndef VECCORE_CUDA
  void Draw(const char *option);
#endif
  bool Resample();

  bool Prune();

VECCORE_ATT_HOST_DEVICE
  int SizeOf() const;
  void Compact();
VECCORE_ATT_HOST_DEVICE
  void RebuildClass();
VECCORE_ATT_HOST_DEVICE
  static int SizeOfStore();
  static int MakeCompactBuffer(char* &b);
VECCORE_ATT_HOST_DEVICE
  static void RebuildStore(char *b);
#ifdef MAGIC_DEBUG
VECCORE_ATT_HOST_DEVICE
  int GetMagic() const {return fMagic;}
#endif
#ifndef VECCORE_CUDA
  static int NLdElems() { return fNLdElems; }

  static TEFstate *Element(int i) {
    if (i < 0 || i >= fNLdElems)
      return 0;
    return fElements[i];
  }
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
 VECCORE_ATT_DEVICE 
 static int NLdElems() { return fEFNLdElemsDev; }
 VECCORE_ATT_DEVICE 
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

#ifdef VECCORE_CUDA
  VECCORE_ATT_HOST_DEVICE
  static TEFstate *GetElement(int z, int a = 0);
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
  VECCORE_ATT_DEVICE TEFstate **GetElements() { return fEFElementsDev; }
#else
  TEFstate **GetElements() { return fEFElementsHost; }
#endif
#else
#ifdef USE_ROOT
  static TEFstate *GetElement(int z, int a = 0, TFile *f = 0);
#endif
  static TEFstate **GetElements() { return fElements; }
#endif
VECCORE_ATT_HOST_DEVICE
  static TPDecay* GetDecayTable() {
#ifndef VECCORE_CUDA
    return fDecay;
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
    return fDecayDev;
#else
    return fDecayHost;
#endif
#endif
  }
VECCORE_ATT_HOST_DEVICE
  static void SetDecayTable(TPDecay *decayTable) {
#ifndef VECCORE_CUDA
    fDecay = decayTable;
#else
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
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

#ifndef VECCORE_CUDA
  static int fNLdElems;              //! number of loaded elements
  static TEFstate *fElements[NELEM]; //! databases of elements
  static TPDecay  *fDecay;           //! decay table
#endif
#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif
VECCORE_ATT_HOST_DEVICE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) fEGrid) % sizeof(fEGrid[0]) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fEGrid misaligned\n");isaligned=false;}
  if(((unsigned long) &fAtcm3) % sizeof(fAtcm3) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fAtcm3 misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmin) % sizeof(fEmin) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fEmin misaligned\n");isaligned=false;}
  if(((unsigned long) &fEmax) % sizeof(fEmax) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fEmax misaligned\n");isaligned=false;}
  if(((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fEilDelta misaligned\n");isaligned=false;}
  if(((unsigned long) &fDens) % sizeof(fDens) != 0) { Geant::Error("TEFstate::CheckAlign","%s","%s","fDens misaligned\n");isaligned=false;}
  if(((unsigned long) &fEle) % sizeof(fEle) != 0) { Geant::Error("TEFstate::CheckAlign","%s","%s","fEle misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fNEbins misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEFstat) % sizeof(fNEFstat) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fNEFstat misaligned\n");isaligned=false;}
  if(((unsigned long) &fNRpart) % sizeof(fNRpart) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fNRpart misaligned\n");isaligned=false;}
  for(auto i=0; i< fNRpart; ++i)
    if(((unsigned long) fPFstateP[i]) % sizeof(double) != 0) { Geant::Error("TEFstate::CheckAlign","fPFstateP[%d misaligned\n",i);isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fMagic misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TEFstate::CheckAlign","%s","fStore misaligned\n");isaligned=false;}
  return isaligned;
}

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
  ClassDefNV(TEFstate, 4) // Element X-secs
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1];   // Pointer to the compact store part of the class
};

#endif
