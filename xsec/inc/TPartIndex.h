// Author: Federico Carminati   27/05/13

/*************************************************************************
 * Copyright (C) 1995-2000, fca                                          *
 * All rights reserved.                                                  *
 *                                                                       *
 *************************************************************************/

#ifndef TPartIndex_H
#define TPartIndex_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TPartIndex                                                           //
//                                                                      //
// Particle index singleton for various particle translation functions  //
//                                                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifdef USE_ROOT
#include "Rtypes.h"
#endif
#ifdef USE_VECGEOM_NAVIGATOR
#include "materials/Particle.h"
#include "base/Map.h"
using vecgeom::Particle;
#else
#include <map>
#endif
#include "Geant/Typedefs.h"

#define DICLEN 12  // Number of process cross sections
#define FNPROC 18  // Number of total processes
#define FNPART 464 // Total number of particles
#define NELEM 118  // Total number of materials

enum GVproc {
  kTransport,
  kMultScatt,
  kIonisation,
  kDecay,
  kinElastic,
  kElastic,
  kRestCapture,
  kBrehms,
  kPairProd,
  kAnnihilation,
  kCoulombScatt,
  kPhotoel,
  kCompton,
  kConversion,
  kCapture,
  kKiller,
  kTotal
};

class TPartIndex {

public:
  static TPartIndex *I() {
     if (!fgPartIndex) {
#ifdef USE_VECGEOM_NAVIGATOR
      Particle::CreateParticles();
#endif
      fgPartIndex = new TPartIndex();
     }
    return fgPartIndex;
  }
  TPartIndex();
  virtual ~TPartIndex();

  static const char *ClassName() { return "TPartIndex"; }

  // Database version
  int Version() const { return fVersion; }
  int VersionMajor() const { return fVersion / 1000 / 1000; }
  int VersionMinor() const { return fVersion / 1000 - VersionMajor() * 1000; }
  int VersionSub() const { return fVersion - VersionMajor() * 1000000 - VersionMinor() * 1000; }

  // Process name <- process index
  const char *ProcName(int proc) const;
  // Process index <- Process name
  int ProcIndex(const char *reac) const {
    int nr = fgNProc;
    while (nr--)
      if (!strcmp(reac, fgPrName[nr]))
        break;
    return nr;
  }

  // Process index <- G4 process*1000+subprocess
  int ProcIndex(int proccode) const;
  // Process index <- G4 process*1000+subprocess
  static int ProcCode(int procindex) /* const */ { return fgPCode[procindex]; }

  static short NProc() /* const */ { return fgNProc; }

  // Fill the particle table
  void SetPartTable(const int *vpdg, int np);

  // PDG code <- GV particle number
  int PDG(int i) const { return fPDG[i]; }
  // PDG code <- particle name
  int PDG(const char *pname) const;
// Particle name <- GV particle number
#ifdef USE_VECGEOM_NAVIGATOR
  const char *PartName(int i) const { return Particle_t::GetParticle(fPDG[i]).Name(); }
#else
  const char *PartName(int i) const { return TDatabasePDG::Instance()->GetParticle(fPDG[i])->GetName(); }
#endif

  // Get the particle from the GeantV code
  const Particle_t *GetParticle(int gvcode) const { return fGVParticle[gvcode]; }

  // GV particle index <- PDG code
  int PartIndex(int pdg) const;

  // GV particle index <- particle name
  int PartIndex(const char *partname) const { return PartIndex(PDG(partname)); }
  // Number of particles
  int NPart() const { return fNPart; }

  // Number of particles with reactions
  void SetNPartReac(int np) { fNpReac = np; }
  void SetNPartCharge(int nc) { fNpCharge = nc; }
  int NPartReac() const { return fNpReac; }
  int NPartCharge() const { return fNpCharge; }
#ifndef USE_VECGEOM_NAVIGATOR
  TDatabasePDG *DBPdg() const { return fDBPdg; }
#endif

  void SetEnergyGrid(double emin, double emax, int nbins);
  int NEbins() const { return fNEbins; }
  double Emin() const { return fEGrid[0]; }
  double Emax() const { return fEGrid[fNEbins - 1]; }
  double EilDelta() const { return fEilDelta; }
  const double *EGrid() const { return fEGrid; }

  static const char *EleSymb(int z) { return fgEleSymbol[z - 1]; }
  static const char *EleName(int z) { return fgEleName[z - 1]; }
  static float WEle(int z) { return fgWElem[z - 1]; }
  static int NElem() { return fgNElem; }

  void Print(const char *option = "") const;
  // approximated formula for nuclear mass computation; for handling fragments
  double GetAprxNuclearMass(int Z, int A);
#ifdef USE_VECGEOM_NAVIGATOR
  void SetPDGToGVMap(vecgeom::map<int, int> &theMap);
#else
  void SetPDGToGVMap(std::map<int, int> &theMap);
#endif
  // only for e-,e+,gamma and proton
  int GetSpecGVIndex(int indx) { return fSpecGVIndices[indx]; }

private:
  TPartIndex(const TPartIndex &);            // Not implemented
  TPartIndex &operator=(const TPartIndex &); // Not implemented

  static TPartIndex *fgPartIndex;

  const int fVersion = 1000002;

  static const int fgNProc = FNPROC;   // Number of processes
  static const char *fgPrName[FNPROC]; // Process name
  static const short fgPCode[FNPROC];  // G4 process codes

  static const int fgNElem = NELEM;      // Number of Elements
  static const char *fgEleSymbol[NELEM]; // Symbol of Element
  static const char *fgEleName[NELEM];   // Name of Element
  static const float fgWElem[NELEM];     // Weight of a mole in grams

  int fNPart;    // Total number of particles
  int *fPDG;     // [fNPart] PDG code of all part
  int fNpReac;   // Number of particles with reactions
  int fNpCharge; // Number of particles with reactions

  int fNEbins;      // number of bins of common energy grid
  double fEilDelta; // Inverse log delta of common energy grid
  double *fEGrid;   // [fNEbins] Common energy grid

#ifndef USE_VECGEOM_NAVIGATOR
  TDatabasePDG *fDBPdg; // Pointer to the augmented pdg database
  std::map<int, int> fPDGToGVMap;              // PDG->GV code map
#else
  vecgeom::map<int, int> fPDGToGVMap;              // PDG->GV code map
#endif
  int fSpecGVIndices[4];                       // store GV codes of e-,e+,gamma and proton
  std::vector<const Particle_t *> fGVParticle; // direct access to particles via GV index

#ifdef USE_ROOT
#ifdef USE_VECGEOM_NAVIGATOR
  ClassDef(TPartIndex, 100) // Particle Index
#else
  ClassDef(TPartIndex, 2) // Particle Index
#endif
#endif
};

#endif
