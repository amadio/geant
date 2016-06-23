//===--- GeantTrack.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTrack.h
 * @brief Implementation of track for GeantV prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK
#define GEANT_TRACK
#define NEW_NAVIGATION

#include "Geant/Config.h"
#include "Geant/Math.h"

#include <climits>
#include <float.h>
#include <atomic>
#include "Geant/Typedefs.h"

#ifndef GEANT_ALIGN_PADDING
#define GEANT_ALIGN_PADDING 64
#endif

#ifndef VECCORE_BITSET_H
#include "base/BitSet.h"
typedef veccore::BitSet BitSet;
#endif

#ifdef GEANT_CUDA
#ifndef GEANT_CUDAUTILS_H
#include "GeantCudaUtils.h"
#endif
#include "backend/cuda/Interface.h"
#endif

/**
 * @enum TrackStatus_t
 */
enum TrackStatus_t { kAlive, kKilled, kInFlight, kBoundary, kExitingSetup, kPhysics, kPostponed, kNew };

/**
 * @enum TransportAction_t
 */
enum TransportAction_t {
  kDone = 0,     /** Return immediately - no tracks left */
  kPostpone = 1, /** return imediately and postpone whatever tracks left */
  kSingle = 2,   /** perform remaining loop in single track mode */
  kVector = 3    /** perform remaining loop in vectorized mode */
};

/**
 * @enum Species_t
 */
enum Species_t { kHadron, kLepton };

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

GEANT_DECLARE_CONSTANT(double, gTolerance);

class GeantTaskData;

/**
 * @brief Class GeantTrack
 */
class GeantTrack {
public:
  int fEvent;            /** Event number */
  int fEvslot;           /** Event slot */
  int fParticle;         /** Index of corresponding particle */
  int fMother;           /** Index of mother particle */
  int fPDG;              /** Particle pdg code */
  int fGVcode;           /** GV particle code */
  int fEindex;           /** Element index */
  int fCharge;           /** Particle charge */
  int fProcess;          /** Current process */
  int fNsteps;           /** Number of steps made */
  Species_t fSpecies;    /** Particle species */
  TrackStatus_t fStatus; /** Track status */
  double fMass;          /** Particle mass */
  double fXpos;          /** X position */
  double fYpos;          /** Y position */
  double fZpos;          /** Z position */
  double fXdir;          /** X direction */
  double fYdir;          /** Y direction */
  double fZdir;          /** Z direction */
  double fP;             /** Momentum */
  double fE;             /** Energy */
  double fTime;          /** Time */
  double fEdep;          /** Energy deposition in the step */
  double fPstep;         /** Selected physical step */
  double fStep;          /** Current step */
  double fSnext;         /** Straight distance to next boundary */
  double fSafety;        /** Safe distance to any boundary */
  double fNintLen;       /** Number of interaction lenghts traveled in last step */
  double fIntLen;        /** Cumulated interaction length since last discrete process */
  bool fBoundary;        /** True if starting from boundary */
  bool fPending;
  VolumePath_t *fPath;
  VolumePath_t *fNextpath;

public:
  /** @brief GeantTrack constructor  */
  GeantTrack();

  /**
   * @brief GeantTrack copy constructor
   */
  GeantTrack(const GeantTrack &other);

  /** @brief Operator = */
  GeantTrack &operator=(const GeantTrack &other);

  /**
  * @brief GeantTrack parametrized constructor
  *
  * @param ipdg ??????
  */
  GeantTrack(int ipdg);

  /**
  * @brief GeantTrack parametrized constructor
  *
  * @param ipdg ??????
  * @param maxdepth ??????
  */
  GEANT_CUDA_BOTH_CODE
  GeantTrack(int ipdg, int maxdepth);

  /** @brief GeantTrack destructor */
  GEANT_CUDA_BOTH_CODE
  ~GeantTrack();

  /** @brief Function that return beta value */
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  double Beta() const { return fP / fE; }

  /** @brief Function that return charge value */
  GEANT_CUDA_BOTH_CODE
  int Charge() const { return fCharge; }

  /** @brief Function that return curvature. To be changed when handling properly field*/
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  double Curvature(double Bz) const {
    // Curvature
    constexpr double kB2C = -0.299792458e-3;
    constexpr double kTiny = 1.E-50;
    double qB = fCharge * Bz;
    if (qB < kTiny) return kTiny;
    return fabs(kB2C * qB / (Pt() + kTiny));
  }

  /** @brief Function that return pointer to X direction value */
  const double *Direction() const { return &fXdir; }

  /** @brief Function that return X direction value */
  double DirX() const { return fXdir; }

  /** @brief Function that return Y direction value */
  double DirY() const { return fYdir; }

  /** @brief Function that return Z direction value */
  double DirZ() const { return fZdir; }

  /** @brief Function that return energy value */
  double E() const { return fE; }

  /** @brief Function that return energy deposition value */
  double Edep() const { return fEdep; }

  /** @brief Function that return event number */
  int Event() const { return fEvent; }

  /** @brief Function that return slot number */
  int EventSlot() const { return fEvslot; }

  /** @brief Function that return true if starting from boundary */
  bool FromBoundary() const { return fBoundary; }

  /** @brief Function that return GV particle code */
  int GVcode() const { return fGVcode; }

  /** @brief Function that return element index */
  int EIndex() const { return fEindex; }

  /** @brief Function that return gamma value*/
  GEANT_CUDA_BOTH_CODE
  double Gamma() const { return fMass ? fE / fMass : DBL_MAX; }

  /** @brief Function that return selected physical step */
  double GetPstep() const { return fPstep; }

  /** @brief Function that return volume */
  Volume_t const*GetVolume() const;

  /** @brief Function that return next volume */
  Volume_t const*GetNextVolume() const;

  /** @brief Function that return material */
  Material_t *GetMaterial() const;

  /** @brief Function that returns the current path (NavigationState) of the track */
  VolumePath_t *GetPath() const { return fPath; }

  /** @brief Function that return next path (NavigationState) of the track */
  VolumePath_t *GetNextPath() const { return fNextpath; }

  /** @brief Function that return number of physical step made */
  int GetNsteps() const { return fNsteps; }

  /** @brief Function that return physical step */
  double GetStep() const { return fStep; }

  /** @brief Function that return time traveled in the step */
  double TimeStep(double step) const { return fE*step/fP; }

  /** @brief Function that return straight distance to next boundary */
  double GetSnext() const { return fSnext; }

  /** @brief Function that return safe distance to any boundary */
  double GetSafety() const { return fSafety; }

  /** @brief Function that return number of interaction lengths travelled since last step */
  double GetNintLen() const { return fNintLen; }

  /** @brief Function that return interaction length since last discrete process */
  double GetIntLen() const { return fIntLen; }

  /** @brief Function that check if track is alive */
  bool IsAlive() const { return (fStatus != kKilled); }

  /** @brief Function that check if track is on boundary */
  bool IsOnBoundary() const { return (fStatus == kBoundary); }

  /** @brief  Check direction normalization within tolerance */
  bool IsNormalized(double tolerance = 1.E-8) const;

  /** @brief Function that set status killed to track */
  void Kill() { fStatus = kKilled; }

  /** @brief Function that return mass value */
  double Mass() const { return fMass; }

  /** @brief Function to normalize direction */
  void Normalize() __attribute__((always_inline)) {
    double norm = 1. / Math::Sqrt(fXdir * fXdir + fYdir * fYdir + fZdir * fZdir);
    fXdir *= norm;
    fYdir *= norm;
    fZdir *= norm;
  }

  /** @brief Function that return momentum value */
  double P() const { return fP; }

  /** @brief Function that return momentum X component */
  double Px() const { return fP * fXdir; }

  /** @brief Function that return momentum Y component */
  double Py() const { return fP * fYdir; }

  /** @brief Function that return momentum Z component */
  double Pz() const { return fP * fZdir; }

  /** @brief Function that return module momentum's value */
  double Pt() const { return fP * Math::Sqrt(fXdir * fXdir + fYdir * fYdir); }

  /** @brief Function that return index of corresponding particle */
  int Particle() const { return fParticle; }

  /** @brief Function that returns index of mother particle */
  int Mother() const { return fMother; }

  /** @brief Function that set status pending to track */
  bool Pending() const { return fPending; }

  /** @brief Function that return particle pdg code */
  int PDG() const { return fPDG; }

  /** @brief Function that return current process */
  int Process() const { return fProcess; }
  const double *Position() const { return &fXpos; }

  /** @brief Function that return X position */
  double PosX() const { return fXpos; }

  /** @brief Function that return Y position */
  double PosY() const { return fYpos; }

  /** @brief Function that return Z position */
  double PosZ() const { return fZpos; }

  /** @brief Print function */
  void Print(const char *location) const;

  /** Function that return particle species */
  Species_t Species() const { return fSpecies; }

  /** Function that return track status */
  TrackStatus_t Status() const { return fStatus; }

  /** Function that return time */
  double Time() const { return fTime; }

  /** Clear function */
  void Clear(const char *option = "");

  /** @brief Function that return X coordinate */
  double X() const { return fXpos; }

  /** @brief Function that return Y coordinate */
  double Y() const { return fYpos; }

  /** @brief Function that return Z coordinate */
  double Z() const { return fZpos; }

  /**
   * @brief Function that set event number
   *
   * @param event Event that should be set as fEvent
   */
  void SetEvent(int event) { fEvent = event; }

  /**
   * @brief Function that set event slot number
   *
   * @param slot Event slot that should be set as fEvslot
   */
  void SetEvslot(int slot) { fEvslot = slot; }

  /**
   * @brief Function that set particle index
   *
   * @param particle Particle that should be set as fParticle
   */
  void SetParticle(int particle) { fParticle = particle; }

  /**
   * @brief Function that sets mother index
   *
   * @param mother Particle that should be set as fMother
   */
  void SetMother(int mother) { fMother = mother; }

  /**
   * @brief Function that set particle pdg code
   *
   * @param pdg Particle pdg code that should be set as fPDG
   */
  void SetPDG(int pdg) { fPDG = pdg; }

  /**
   * @brief Function that set GV particle code
   *
   * @param gVcode GV particle code that should be set as fGVcode
   */
  void SetGVcode(int gVcode) { fGVcode = gVcode; }

  /**
   * @brief Function that set element index
   *
   * @param ind Element index that should be set as fEindex
   */
  void SetEindex(int ind) { fEindex = ind; }

  /**
   * @brief Function that set charge
   *
   * @param charge Charge that should be set as fCharge
   */
  void SetCharge(int charge) { fCharge = charge; }

  /**
   * @brief Function that set process
   *
   * @param process Process that should be set as fProcess
   */
  void SetProcess(int process) { fProcess = process; }

  /**
   * @brief Function that set current step
   *
   * @param nsteps Current step hat should be set as fNsteps
   */
  void SetNsteps(int nsteps) { fNsteps = nsteps; }

  /**
   * @brief Function that set current species
   *
   * @param species Current species hat should be set as fSpecies
   */
  void SetSpecies(Species_t species) { fSpecies = species; }

  /**
   * @brief Function that set track status
   *
   * @param status Current track status that should be set as fStatus
   */
  void SetStatus(TrackStatus_t &status) { fStatus = status; }

  /**
   * @brief Function that set mass
   *
   * @param mass Current mass that should be set as fMass
   */
  void SetMass(double mass) { fMass = mass; }

  /**
   * @brief Function that set X, Y, Z positions
   *
   * @param x X position
   * @param y Y position
   * @param z Z position
   */
  void SetPosition(double x, double y, double z) {
    fXpos = x;
    fYpos = y;
    fZpos = z;
  }

  /**
   * @brief [Function that set X, Y, Z directions
   *
   * @param dx X direction
   * @param dy Y direction
   * @param dz Z direction
   */
  void SetDirection(double dx, double dy, double dz) {
    fXdir = dx;
    fYdir = dy;
    fZdir = dz;
  }

  /**
   * @brief Function that set momentum
   *
   * @param p Current momentum should be set as fP
   */
  void SetP(double p) { fP = p; }

  /**
   * @brief Function that set energy
   *
   * @param e Current E should be set as fE
   */
  void SetE(double e) { fE = e; }

  /**
   * @brief Function that set time
   *
   * @param time Current time should be set as fTime
   */
  void SetTime(double time) { fTime = time; }

  /**
   * @brief Function that set energy deposition
   *
   * @param edep Current energy deposition should be set as fEdep
   */
  void SetEdep(double edep) { fEdep = edep; }

  /**
   * @brief Function that set current physical step
   *
   * @param pstep Current physical step should be set as fPstep
   */
  void SetPstep(double pstep) { fPstep = pstep; }

  /**
   * @brief Function that set straight distance to next boundary
   *
   * @param snext Straight distance to next boundary should be set as fSnext
   */
  void SetSnext(double snext) { fSnext = snext; }

  /**
   * @brief Function that set safe distance to any boundary
   *
   * @param safety Safe distance to any boundary hould be set as fSafety
   */
  void SetSafety(double safety) { fSafety = safety; }

  /**
   * @brief Function that set number of interaction lengths traveled
   *
   * @param nradlen Number of interaction lengths travelled
   */
  void SetNintLen(double nradlen) { fNintLen = nradlen; }

  /**
   * @brief Function that set interaction length since last discrete process
   *
   * @param intlen Interaction length
   */
  void SetIntLen(double intlen) { fIntLen = intlen; }

  /**
   * @brief Function that set starting from boundary flag
   *
   * @param flag Flag that is true if starting from boundary
   */
  void SetFrombdr(bool flag) { fBoundary = flag; }

  /**
   * @brief Function that set pending status
   *
   * @param flag Flag that should be set pending
   */
  void SetPending(bool flag) { fPending = flag; }

  /**
   * @brief Function that set next path
   *
   * @param path Volume path
   */
  void SetPath(VolumePath_t const *const path);

  /**
   * @brief Function that set next volume path
   *
   * @param path Volume path
   */
  void SetNextPath(VolumePath_t const *const path);

  //  ClassDefNV(GeantTrack, 1) // The track
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
