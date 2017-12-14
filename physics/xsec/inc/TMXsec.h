// Author: Federico Carminati   27/05/13
/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/
#ifndef TMXsec_H
#define TMXsec_H
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TMXSec                                                               //
//                                                                      //
// X-section for GV per material                                        //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#ifndef VECCORE_CUDA
#include <vector>
#else 
#include "base/Vector.h"
#include "base/Map.h"
using vecgeom::Vector;
using vecgeom::pair;
#endif

#include "Geant/Config.h"
#include "GeantTrack.h"

#include <TEXsec.h>

class TPDecay;
#include "GeantFwd.h"

class TMXsec {
public:
  using GeantTrack = Geant::GeantTrack;
  using GeantTrack_v = Geant::GeantTrack_v;
  using GeantTaskData = Geant::GeantTaskData;
  using TrackVec_t = Geant::TrackVec_t;

  TMXsec();
  TMXsec(const char *name, const char *title, const int z[], const int a[], const float w[], int nel, float dens,
         bool weight = false, const TPDecay *decaytable = 0);
  virtual ~TMXsec();
  static const char *ClassName() { return "TMXsec"; }
  const char *GetName() const { return fName; }
  const char *GetTitle() const { return fTitle; }
  float Xlength(int part, float en, double ptot);
  //   bool Xlength_v(int npart, const int part[], const float en[], double lam[]);
  float DEdx(int part, float en);
  //   bool DEdx_v(int npart, const int part[], const float en[], float de[]);
  VECCORE_ATT_HOST_DEVICE
  float Range(int part, float en);
  VECCORE_ATT_HOST_DEVICE
  double InvRange(int part, float step);

  VECCORE_ATT_HOST_DEVICE
  void Eloss(int ntracks, GeantTrack_v &tracks,GeantTaskData *td);
  VECCORE_ATT_HOST_DEVICE
  void ElossSingle(int itrack, GeantTrack_v &tracks,GeantTaskData *td);
  void ProposeStep(int ntracks, GeantTrack_v &tracks, GeantTaskData *td);
  void ProposeStepSingle(int itr, GeantTrack_v &tracks, GeantTaskData *td);
  void SampleInt(int ntracks, GeantTrack_v &tracksin, GeantTaskData *td);
  void SampleSingleInt(int itr, GeantTrack_v &tracksin, GeantTaskData *td);

//=== N E W   I N T E R F A C E S ===//
  VECCORE_ATT_HOST_DEVICE
  void Eloss(GeantTrack &track, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void Eloss(TrackVec_t &tracks, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void ProposeStep(GeantTrack &track, GeantTaskData *td);
  
  VECCORE_ATT_HOST_DEVICE
  void ProposeStep(TrackVec_t &tracks, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void SampleInt(GeantTrack &track, GeantTaskData *td);

  VECCORE_ATT_HOST_DEVICE
  void SampleInt(TrackVec_t &tracks, GeantTaskData *td);
//===================================//

  VECCORE_ATT_HOST_DEVICE
  float MS(int ipart, float energy);

  TEXsec *SampleInt(int part, double en, int &reac, double ptotal);
  int SampleElement(GeantTaskData *td); // based on # atoms/vol. for the prototype
  int SampleElement();                  // based on # atoms/vol. for Geant4 with tab.phys.

  int SelectElement(int pindex, int rindex, double energy);

  //   static bool Prune();
  void Print(const char *opt = "") const;

private:
  TMXsec(const TMXsec &);            // Not implemented
  TMXsec &operator=(const TMXsec &); // Not implemented

  char fName[32];       // cross section name
  char fTitle[128];     // cross section title
  int fNEbins;          // number of energy bins
  int fNTotXL;          // dimension of fTotXL
  int fNCharge;         // dimension of tables for charged particles
  int fNRelXS;          // dimension of fRelXS
  double fEilDelta;     // logarithmic energy delta
  const double *fEGrid; //! Energy grid

  int fNElems;                                            // Number of elements
  TEXsec **fElems;                                        // [fNElems] List of elements composing this material
  float *fTotXL;                                          // [fNTotXL] Total x-sec for this material
  float *fRelXS;                                          // [fNRelXS] Relative x-sec for this material
  float *fDEdx;                                           // [fNCharge] Ionisation energy loss for this material
  float *fMSangle;                                        // [fNCharge] table of MS average angle
  float *fMSansig;                                        // [fNCharge] table of MS sigma angle
  float *fMSlength;                                       // [fNCharge] table of MS average lenght correction
  float *fMSlensig;                                       // [fNCharge] table of MS sigma lenght correction
  double *fRatios;                                        // [fNElems]  relative #atoms/volume; normalized
  float *fRange;                                          // [fNCharge] ranges of the particle in this material
  const TPDecay *fDecayTable;                             // pointer to the decay table
#ifndef VECCORE_CUDA
  std::vector<std::pair<float, double>> **fInvRangeTable; // [fNCharge]
#else
  Vector<vecgeom::pair<float, double>> **fInvRangeTable; // [fNCharge]

#endif
};

#ifdef USE_ROOT
class TOMXsec : public TObject {
public:
  TOMXsec(TMXsec *mxsec) : fMXsec(mxsec) {}
  ~TOMXsec() { delete fMXsec; }
  TMXsec *MXsec() const { return fMXsec; }

private:
  TOMXsec(const TOMXsec &);            // Not implemented
  TOMXsec &operator=(const TOMXsec &); // Not implemented

  TMXsec *fMXsec;
};
#endif
#endif
