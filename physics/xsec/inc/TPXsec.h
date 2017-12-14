// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TPXsec_H
#define TPXsec_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPXSec                                                               //
//                                                                      //
// X-section for GV per particle                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TPartIndex.h"
#include "Geant/Error.h"

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#endif
 
class TPXsec {
public:
  VECCORE_ATT_HOST_DEVICE
  TPXsec();
  VECCORE_ATT_HOST_DEVICE
  TPXsec(int pdg, int nxsec);
  VECCORE_ATT_HOST_DEVICE
  TPXsec(const TPXsec &other);
  VECCORE_ATT_HOST_DEVICE
  virtual ~TPXsec();
  void Print(const char *opt = "") const;
  VECCORE_ATT_HOST_DEVICE
  const char *Name() const { return TPartIndex::I()->PartName(fPDG); }
  bool SetPart(int pdg, int nxsec);
  bool SetPartXS(const float xsec[], const int dict[]);
  bool SetPartIon(const float dedx[]);
  bool SetPartMS(const float angle[], const float ansig[], const float length[], const float lensig[]);
  int PDG() const { return fPDG; }
  VECCORE_ATT_HOST_DEVICE
  float XS(int rindex, double en, bool verbose=false) const;
  bool XS_v(int npart, int rindex, const double en[], double lam[]) const;
  float DEdx(double en) const;
  bool MS(double en, float &ang, float &asig, float &len, float &lsig) const;
  VECCORE_ATT_HOST_DEVICE
  bool Resample();
  bool Prune();
  int SampleReac(double en) const;
  int SampleReac(double en, double randn) const;

  void Dump() const;
  VECCORE_ATT_HOST_DEVICE
  void Interp(double egrid[], float value[], int nbins, double eildelta, int stride, double en, float result[]);

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }
  VECCORE_ATT_HOST_DEVICE
  int SizeOf() const;
  VECCORE_ATT_HOST_DEVICE
  void Compact();
  VECCORE_ATT_HOST_DEVICE
  void RebuildClass();
#ifdef MAGIC_DEBUG
  VECCORE_ATT_HOST_DEVICE
  int GetMagic() const {return fMagic;}
#endif

  VECCORE_ATT_HOST_DEVICE
bool CheckAlign() {
  bool isaligned=true;
  if(((unsigned long) &fPDG) % sizeof(fPDG) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fPDG misaligned\n");isaligned=false;}
  if(((unsigned long) &fNEbins) % sizeof(fNEbins) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fNEbins misaligned\n");isaligned=false;}
  if(((unsigned long) &fNCbins) % sizeof(fNCbins) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fNCbins misaligned\n");isaligned=false;}
  if(((unsigned long) &fNXsec) % sizeof(fNXsec) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fNXsec misaligned\n");isaligned=false;}
  if(((unsigned long) &fNTotXs) % sizeof(fNTotXs) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fNTotXs misaligned\n");isaligned=false;}
  if(((unsigned long) &fNXSecs) % sizeof(fNXSecs) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fNXSecs misaligned\n");isaligned=false;}
  if(((unsigned long) fEGrid) % sizeof(fEGrid[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fEGrid misaligned\n");isaligned=false;}
  if(((unsigned long) fMSangle) % sizeof(fMSangle[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fMSangle misaligned\n");isaligned=false;}
  if(((unsigned long) fMSansig) % sizeof(fMSansig[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fMSansig misaligned\n");isaligned=false;}
  if(((unsigned long) fMSlength) % sizeof(fMSlength[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fMSlength misaligned\n");isaligned=false;}
  if(((unsigned long) fMSlensig) % sizeof(fMSlensig[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fMSlensig misaligned\n");isaligned=false;}
  if(((unsigned long) fdEdx) % sizeof(fdEdx[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fdEdx misaligned\n");isaligned=false;}
  if(((unsigned long) fTotXs) % sizeof(fTotXs[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fTotXs misaligned\n");isaligned=false;}
  if(((unsigned long) fXSecs) % sizeof(fXSecs[0]) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fXSecs misaligned\n");isaligned=false;}
  if(int delta = ((unsigned long) &fEmin) % sizeof(fEmin) != 0) { Geant::Error("TPXsec::CheckAlign","fEmin misaligned %d \n",delta);isaligned=false;}
  if(int delta = ((unsigned long) &fEmax) % sizeof(fEmax) != 0) { Geant::Error("TPXsec::CheckAlign","fEmax misaligned %d \n",delta);isaligned=false;}
  if(int delta = ((unsigned long) &fEilDelta) % sizeof(fEilDelta) != 0) { Geant::Error("TPXsec::CheckAlign","fEilDelta misaligned %d \n",delta);isaligned=false;}
  if(((unsigned long) &fRdict) % sizeof(int) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fRdict misaligned\n");isaligned=false;}
  if(((unsigned long) &fRmap) % sizeof(int) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fRmap misaligned\n");isaligned=false;}
#ifdef MAGIC_DEBUG
  if(((unsigned long) &fMagic) % sizeof(fMagic) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fMagic misaligned\n");isaligned=false;}
#endif
  if(((unsigned long) &fStore) % sizeof(double) != 0) { Geant::Error("TPXsec::CheckAlign","%s","fStore misaligned\n");isaligned=false;}
  return isaligned;
}
#ifdef VECCORE_CUDA
VECCORE_ATT_HOST_DEVICE
char *strncpy(char *dest, const char *src, size_t n)
{
    char *ret = dest;
    do {
        if (!n--)
            return ret;
    } while (*dest++ = *src++);
    while (n--)
        *dest++ = 0;
    return ret;
};
#endif

private:
  TPXsec &operator=(const TPXsec &); // Not implemented

  static int fVerbose; // Controls verbosity level

  int fPDG;             // particle pdg code
  int fNEbins;          // number of energy bins
  int fNCbins;          // number of energy bins for dEdx and MS
  int fNXsec;           // number of reactions
  int fNTotXs;          // tot size of fTotXs
  int fNXSecs;          // tot size of fXSecs
  const double *fEGrid; //![fNEbins] energy grid
  float *fMSangle;      // [fNCbins] table of MS average angle
  float *fMSansig;      // [fNCbins] table of MS sigma angle
  float *fMSlength;     // [fNCbins] table of MS average lenght correction
  float *fMSlensig;     // [fNCbins] table of MS sigma lenght correction
  float *fdEdx;         // [fNCbins] table of dE/dx
  float *fTotXs;        // [fNTotXs] table of total x-sec
  float *fXSecs;        // [fNXSecs] table of partial x-sec
  double fEmin;         // Min energy of the energy grid
  double fEmax;         // Max energy of the energy grid
  double fEilDelta;     // logarithmic energy delta
  int fRdict[FNPROC];   // reaction dictionary from reaction number to position
                        // in the X-sec array
  int fRmap[FNPROC];    // reaction map, from reaction position in the X-sec
                        // array to the raction number
#ifdef MAGIC_DEBUG
  const int fMagic = -777777;
#endif

#ifndef VECCORE_CUDA
#ifdef USE_ROOT
  ClassDefNV(TPXsec, 4) // Particle X-secs
#endif
#endif

private:
  alignas(sizeof(double)) char fStore[1];        //! Pointer to the compact data of the class
};

#endif
