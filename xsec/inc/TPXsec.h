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
#ifdef USE_ROOT
#include "Rtypes.h"
#endif

class TPXsec {
public:
  TPXsec();
  TPXsec(int pdg, int nxsec);
  TPXsec(const TPXsec &other);
  virtual ~TPXsec();
  void Print(const char *opt = "") const;
  const char *Name() const { return TPartIndex::I()->PartName(fPDG); }
  bool SetPart(int pdg, int nxsec);
  bool SetPartXS(const float xsec[], const int dict[]);
  bool SetPartIon(const float dedx[]);
  bool SetPartMS(const float angle[], const float ansig[], const float length[], const float lensig[]);
  int PDG() const { return fPDG; }
  float XS(int rindex, double en) const;
  bool XS_v(int npart, int rindex, const double en[], double lam[]) const;
  float DEdx(double en) const;
  bool MS(double en, float &ang, float &asig, float &len, float &lsig) const;
  bool Resample();
  bool Prune();
  int SampleReac(double en) const;
  int SampleReac(double en, double randn) const;

  void Dump() const;
  void Interp(double egrid[], float value[], int nbins, double eildelta, int stride, double en, float result[]);

  static void SetVerbose(int verbose) { fVerbose = verbose; }
  static int GetVerbose() { return fVerbose; }
  int SizeOf() const;
  void Compact();
  void RebuildClass();
#ifdef MAGIC_DEBUG
  int GetMagic() const {return fMagic;}
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
#ifdef USE_ROOT
  ClassDefNV(TPXsec, 3) // Particle X-secs
#endif
private:
  float *fStore;        //! Pointer to the compact data of the class
};

#endif
