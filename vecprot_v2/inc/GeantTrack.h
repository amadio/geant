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

#include "base/Vector3D.h"
#include "Geant/Config.h"
#include "Geant/Math.h"

#include <functional>
#include <climits>
#include <float.h>
#include <atomic>
#include <mutex>
#include "Geant/Typedefs.h"

#ifndef GEANT_ALIGN_PADDING
#define GEANT_ALIGN_PADDING 64
#endif

#ifndef VECCORE_BITSET_H
#include "base/BitSet.h"
typedef veccore::BitSet BitSet;
#endif

#ifndef VECCORE_CUDA
#ifdef GEANT_USE_NUMA
#include "NumaAllocator.h"
#endif
#endif

#ifdef GEANT_CUDA_ENABLED
#include "GeantCudaUtils.h"
#include "backend/cuda/Interface.h"
#endif

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

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

/** Basket simulation stages. */
enum ESimulationStage {
#ifdef USE_REAL_PHYSICS
  kPreStepStage         = 0, // Actions at the beginning of the step
  kComputeIntLStage,         // Physics interaction length computation stage
  kGeometryStepStage,        // Compute geometry transport length
  kPrePropagationStage,      // Special msc stage for step limit phase
//  kGeometryStepStage,        // Compute geometry transport length
  kPropagationStage,         // Propagation in field stage
  kPostPropagationStage,     // Special msc stage for along-step action stage
//  kMSCStage,               // Multiple scattering stage
  kAlongStepActionStage,     // Along step action stage (continuous part of the inetraction)
  kPostStepActionStage,      // Post step action stage (discrete part of the inetraction)
// kAtRestActionStage,       // At-rest action stage (at-rest part of the inetraction)
  kSteppingActionsStage      // User actions
#else
  kPreStepStage         = 0, // Actions at the beginning of the step
  kXSecSamplingStage,        // Propose physics step by sampling total Xsec
  kGeometryStepStage,        // Compute geometry transport length
  kPropagationStage,         // Propagation in field stage
  //  kMSCStage,                 // Multiple scattering stage
  kContinuousProcStage,      // Continuous processes stage
  kDiscreteProcStage,        // Discrete processes stage
  kSteppingActionsStage      // User actions
#endif
};

class GeantTaskData;
class GeantTrack;
class TrackDataMgr;

using InplaceConstructFunc_t = std::function<void(void*)>;
using PrintUserDataFunc_t = std::function<void(void*)>;
using Vector3 = vecgeom::Vector3D<double>;

struct StepPointInfo {
  Vector3 fPos;          /** Position */
  Vector3 fDir;          /** Direction */
  double fE;             /** Energy */
  double fP;             /** Momentum */
  double fTime;          /** Time */
  VolumePath_t *fPath;   /** Paths for the particle in the geometry */
};

GEANT_DECLARE_CONSTANT(double, gTolerance);

#ifndef VECCORE_CUDA
#ifdef GEANT_USE_NUMA
typedef NumaAllocator<GeantTrack*> TrackAllocator_t;
typedef std::vector<GeantTrack *, TrackAllocator_t> TrackVec_t;
#else
typedef vector_t<GeantTrack *> TrackVec_t;
#endif
#else
typedef vecgeom::Vector<GeantTrack *> TrackVec_t;
#endif

/**
 * @brief Class GeantTrack
 */
class GeantTrack {
public:
  int fEvent = -1;           /** Event number */
  int fEvslot = -1;          /** Event slot */
  int fParticle = -1;        /** Index of corresponding particle */
  int fPrimaryIndx = -1;     /** Index of the primary particle in the cuurent event */
  int fMother = 0;           /** Index of mother particle */
  int fPDG = 0;              /** Particle pdg code */
  int fGVcode = 0;           /** GV particle code */
  int fEindex = 0;           /** Element index */
  int fBindex = 0;           /** Index in the track block */
  int fCharge = 0;           /** Particle charge */
  int fProcess = -1;         /** Current process */
  int fNsteps = 0;           /** Number of steps made */
  int fMaxDepth = 0;         /** Maximum geometry depth */
  int fStage = 0;            /** Simulation stage */
  int fGeneration = 0;       /** Track generation: 0=primary */
  Species_t fSpecies = kHadron;   /** Particle species */
  TrackStatus_t fStatus = kAlive; /** Track status */
  double fMass = 0;          /** Particle mass */
  double fXpos = 0;          /** X position */
  double fYpos = 0;          /** Y position */
  double fZpos = 0;          /** Z position */
  double fXdir = 0;          /** X direction */
  double fYdir = 0;          /** Y direction */
  double fZdir = 0;          /** Z direction */
  double fP = 0;             /** Momentum */
  double fE = 0;             /** Energy */
  double fTime = 0;          /** Time */
  double fEdep = 0;          /** Energy deposition in the step */
  double fPstep = 1.e+20;    /** Selected physical step */
  double fStep = 0;          /** Current step */
  double fSnext = 0;         /** Straight distance to next boundary */
  double fSafety = 0;        /** Safe distance to any boundary */
  double fNintLen = 0;       /** Number of interaction lenghts traveled in last step */
  double fIntLen = 0;        /** Cumulated interaction length since last discrete process */
  bool fBoundary = false;    /** True if starting from boundary */
  bool fPending = false;     /** Track pending to be processed */
  bool fOwnPath = false;     /** Marker for path ownership */
  bool fIsOnBoundaryPreStp = false;   // to indicate that the particle was on boundary at the pre-step pint
  Volume_t const *fVolume = nullptr; /** Current volume the particle is in */

  // max number of physics processesper particle is assumed to be 10!
  static constexpr size_t fNumPhysicsProcess = 10;
  size_t  fPhysicsProcessIndex = -1;  // selected physics process
  double  fPhysicsNumOfInteractLengthLeft[fNumPhysicsProcess];
  double  fPhysicsInteractLength[fNumPhysicsProcess]; // mfp

private:
#ifdef GEANT_NEW_DATA_FORMAT
  StepPointInfo fPreStep;    /** Pre-step state */
  StepPointInfo fPostStep;   /** Post-step state */
#else
  VolumePath_t *fPath = nullptr;     /** Paths for the particle in the geometry */
  VolumePath_t *fNextpath = nullptr; /** Path for next volume */
#endif

public:
char *fExtraData;                   /** Start of user data at first aligned address. This needs to point to an aligned address. */
  // See: implementation of SizeOfInstance. Offsets of user data blocks are with respect to the fExtraData address. We have to pinpoint this
  // because the start address of the block of user data needs to be aligned, which make the offsets track-by-track dependent if they are
  // relative to the start address of the track (which has to be un-aligned for space considerations)
  // The track memory layout will be as follows (* mark alignment padding start, _ mark padding empty space).
  // Sizes are managed by TrackDataMgr singleton.

  // *__|addr______________*_________________*_________________*_________________*_________________*_________________*_________________*
  //  __|++++track_data+++++++__padding______|++user_data_1++__|++user_data_2++_..._____________   |++fPath_pre_step+++|_______________|++fPath_post_step
  //Addresses:
  //    | addr (track start address)
  //                                         | fExtraData = round_up_align(addr + sizeof(GeantTrack))
  //                                                           | round_up_align(fExtraData + sizeof(user_data_1))
  //                                                                                               | round_up_align(fExtraData + TrackDataMgr::fDataOffset))


  /**
  * @brief GeantTrack in place constructor
  *
  * @param addr Start address
  */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(void *addr);

public:
  /** @brief GeantTrack destructor */
  VECCORE_ATT_HOST_DEVICE
  ~GeantTrack() = delete;

  /**
   * @brief GeantTrack copy constructor
   */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(const GeantTrack &other) = delete;

  /** @brief Operator = */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack &operator=(const GeantTrack &other);

  /** @brief Function that return beta value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Beta() const { return fP / fE; }

  /** @brief Function that return charge value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Charge() const { return fCharge; }

  /** @brief Function that return curvature. To be changed when handling properly field*/
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Curvature(double Bz) const {
    // Curvature
    constexpr double kB2C = -0.299792458e-3;
    constexpr double kTiny = 1.E-50;
    double qB = fCharge * Bz;
    if (fabs(qB) < kTiny) return kTiny;
    return fabs(kB2C * qB / (Pt() + kTiny));
  }

  /** @brief Function that return pointer to X direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  const double *Direction() const { return &fXdir; }

  /** @brief Function that return X direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double DirX() const { return fXdir; }

  /** @brief Function that return Y direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double DirY() const { return fYdir; }

  /** @brief Function that return Z direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double DirZ() const { return fZdir; }

  /** @brief Function that return energy value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double E() const { return fE; }

  /** @brief Function that return kinetic energy value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double T() const { return (fE - fMass); }

  /** @brief Function that stops the track depositing its kinetic energy */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Stop() { fEdep = T(); fE = fMass; fP = 0; }

  /** @brief Function that return energy deposition value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Edep() const { return fEdep; }

  /** @brief Function that return event number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Event() const { return fEvent; }

  /** @brief Function that return slot number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int EventSlot() const { return fEvslot; }

  /** @brief Function that return true if starting from boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool FromBoundary() const { return fBoundary; }

  /** @brief Function that return GV particle code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GVcode() const { return fGVcode; }

  /** @brief Function that return element index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int EIndex() const { return fEindex; }

  /** @brief Function that return index in the track block */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int BIndex() const { return fBindex; }

  /** @brief Function that return gamma value*/
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Gamma() const { return fMass ? fE / fMass : DBL_MAX; }

  /** @brief Function that return selected physical step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetPstep() const { return fPstep; }

  /** @brief Function that return volume */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Volume_t const *GetVolume() const { return fVolume; }

  /** @brief Function that returns current path */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  VolumePath_t *Path() const { return fPath; }

  /** @brief Function that returns next path */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  VolumePath_t *NextPath() const { return fNextpath; }

  /** @brief Function that swaps path and next path */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdateSwapPath() { 
    VolumePath_t *tmp = fNextpath;
    fNextpath = fPath;
    fPath = tmp;
    UpdateVolume();
  }

  /** @brief Function that makes next path have the same content as the current. */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdateSameNextPath() { *fNextpath = *fPath; }
  
  /** @brief Function that updates the current volume the particle is in */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdateVolume() {   
#ifdef USE_VECGEOM_NAVIGATOR
    fVolume = fPath->Top()->GetLogicalVolume();
#else
    fVolume = fPath->GetCurrentNode()->GetVolume();
#endif
  }

  /** @brief Function that return next volume */
  VECCORE_ATT_HOST_DEVICE
  Volume_t const *GetNextVolume() const {
#ifdef USE_VECGEOM_NAVIGATOR
    return ( fNextpath->Top()->GetLogicalVolume() );
#else
    return ( fNextpath->GetCurrentNode()->GetVolume() );
#endif
  }

  /** @brief Function that return material */
  Material_t *GetMaterial() const;

  /** @brief Function that return number of physical step made */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNsteps() const { return fNsteps; }

  /** @brief Function that return physical step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetStep() const { return fStep; }

  /** @brief Function that return time traveled in the step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double TimeStep(double step) const { return fE*step/fP; }

  /** @brief Function that return straight distance to next boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetSnext() const { return fSnext; }

  /** @brief Function that return safe distance to any boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetSafety() const { return fSafety; }

  /** @brief Function that return number of interaction lengths travelled since last step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetNintLen() const { return fNintLen; }

  /** @brief Function that return interaction length since last discrete process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetIntLen() const { return fIntLen; }

  /** @brief Function that check if track is alive */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsAlive() const { return (fStatus != kKilled); }

  /** @brief Function that check if track is on boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsOnBoundary() const { return (fStatus == kBoundary); }

  /** @brief  Check direction normalization within tolerance */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsNormalized(double tolerance = 1.E-8) const {
    double norm = fXdir * fXdir + fYdir * fYdir + fZdir * fZdir;
    if (fabs(1. - norm) > tolerance)
      return false;
    return true;
  }

  /** @brief Function that set status killed to track */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Kill() { fStatus = kKilled; }

  /** @brief Function that return mass value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Mass() const { return fMass; }

  /** @brief Function to normalize direction */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Normalize() {
    double norm = 1. / Math::Sqrt(fXdir * fXdir + fYdir * fYdir + fZdir * fZdir);
    fXdir *= norm;
    fYdir *= norm;
    fZdir *= norm;
  }

  /** @brief Function to make a step along the current direction */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void MakeStep(double step) {
    fPstep -= step;
    fStep += step;
    fSafety -= step;
    fSnext -= step;
    fXpos += step * fXdir;
    fYpos += step * fYdir;
    fZpos += step * fZdir;
  }


  /** @brief Function that return momentum value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double P() const { return fP; }

  /** @brief Function that return momentum X component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Px() const { return fP * fXdir; }

  /** @brief Function that return momentum Y component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Py() const { return fP * fYdir; }

  /** @brief Function that return momentum Z component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Pz() const { return fP * fZdir; }

  /** @brief Function that return module momentum's value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Pt() const { return fP * Math::Sqrt(fXdir * fXdir + fYdir * fYdir); }

  /** @brief Function that return index of corresponding particle */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Particle() const { return fParticle; }

  /** @brief Function that return index of the primary particle in the current event */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int PrimaryParticleIndex() const { return fPrimaryIndx; }

  /** @brief Function that returns index of mother particle */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Mother() const { return fMother; }

  /** @brief Function that set status pending to track */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool Pending() const { return fPending; }

  /** @brief Function that return particle pdg code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int PDG() const { return fPDG; }

  /** @brief Function that return current process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Process() const { return fProcess; }
  const double *Position() const { return &fXpos; }

  /** @brief Function that return X position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double PosX() const { return fXpos; }

  /** @brief Function that return Y position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double PosY() const { return fYpos; }

  /** @brief Function that return Z position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double PosZ() const { return fZpos; }

  /** @brief Print function */
  void Print(const char *msg = "") const;

  /** @brief Print function for a container of tracks */
  static void PrintTracks(TrackVec_t &tracks);

  /** Function that return particle species */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Species_t Species() const { return fSpecies; }

  /** Function that return track status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  TrackStatus_t Status() const { return fStatus; }

  /** Function that return time */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Time() const { return fTime; }

  /** Clear function */
  VECCORE_ATT_HOST_DEVICE
  void Clear(const char *option = "");

  /** Fast reset function */
  VECCORE_ATT_HOST_DEVICE
  void Reset(GeantTrack const &blueprint);

  /** @brief Function that return X coordinate */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double X() const { return fXpos; }

  /** @brief Function that return Y coordinate */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Y() const { return fYpos; }

  /** @brief Function that return Z coordinate */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Z() const { return fZpos; }

  /**
   * @brief Function that set event number
   *
   * @param event Event that should be set as fEvent
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEvent(int event) { fEvent = event; }

  /**
   * @brief Function that set event slot number
   *
   * @param slot Event slot that should be set as fEvslot
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEvslot(int slot) { fEvslot = slot; }

  /**
   * @brief Function that sets the primary particle index in the current event
   *
   * @param primaryindx Index of the primary particle in the current event that this track belongs to
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPrimaryParticleIndex(int primaryindx) { fPrimaryIndx = primaryindx; }

  /** @brief Setter for stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStage(ESimulationStage stage) { fStage = stage; }

  /**
   * @brief Function that sets mother index
   *
   * @param mother Particle that should be set as fMother
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetMother(int mother) { fMother = mother; }

  /**
   * @brief Function that set particle pdg code
   *
   * @param pdg Particle pdg code that should be set as fPDG
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPDG(int pdg) { fPDG = pdg; }

  /**
   * @brief Function that set GV particle code
   *
   * @param gVcode GV particle code that should be set as fGVcode
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetGVcode(int gVcode) { fGVcode = gVcode; }

  /**
   * @brief Function that set element index
   *
   * @param ind Element index that should be set as fEindex
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEindex(int ind) { fEindex = ind; }

  /**
   * @brief Function that set index of the track block
   *
   * @param ind Index to be set
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBindex(int ind) { fBindex = ind; }

  /**
   * @brief Function that set charge
   *
   * @param charge Charge that should be set as fCharge
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetCharge(int charge) { fCharge = charge; }

  /**
   * @brief Function that set process
   *
   * @param process Process that should be set as fProcess
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetProcess(int process) { fProcess = process; }

  /**
   * @brief Function that set current step
   *
   * @param nsteps Current step hat should be set as fNsteps
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNsteps(int nsteps) { fNsteps = nsteps; }

  /**
   * @brief Function that set current species
   *
   * @param species Current species hat should be set as fSpecies
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSpecies(Species_t species) { fSpecies = species; }

  /**
   * @brief Function that set track status
   *
   * @param status Current track status that should be set as fStatus
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStatus(TrackStatus_t &status) { fStatus = status; }

  /**
   * @brief Function that set mass
   *
   * @param mass Current mass that should be set as fMass
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetMass(double mass) { fMass = mass; }

  /**
   * @brief Function that set X, Y, Z positions
   *
   * @param x X position
   * @param y Y position
   * @param z Z position
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
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
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
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
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetP(double p) { fP = p; }

  /**
   * @brief Function that set energy
   *
   * @param e Current E should be set as fE
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetE(double e) { fE = e; }

  /**
   * @brief Function that set time
   *
   * @param time Current time should be set as fTime
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetTime(double time) { fTime = time; }

  /**
   * @brief Function that set energy deposition
   *
   * @param edep Current energy deposition should be set as fEdep
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEdep(double edep) { fEdep = edep; }

  /**
   * @brief Function that set current physical step
   *
   * @param pstep Current physical step should be set as fPstep
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPstep(double pstep) { fPstep = pstep; }

  /**
   * @brief Function that set straight distance to next boundary
   *
   * @param snext Straight distance to next boundary should be set as fSnext
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSnext(double snext) { fSnext = snext; }

  /**
   * @brief Function that set safe distance to any boundary
   *
   * @param safety Safe distance to any boundary hould be set as fSafety
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSafety(double safety) { fSafety = safety; }

  /**
   * @brief Function that set number of interaction lengths traveled
   *
   * @param nradlen Number of interaction lengths travelled
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNintLen(double nradlen) { fNintLen = nradlen; }

  /**
   * @brief Function that set interaction length since last discrete process
   *
   * @param intlen Interaction length
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetIntLen(double intlen) { fIntLen = intlen; }

  /**
   * @brief Function that set starting from boundary flag
   *
   * @param flag Flag that is true if starting from boundary
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetFrombdr(bool flag) { fBoundary = flag; }

  /**
   * @brief Function that set pending status
   *
   * @param flag Flag that should be set pending
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPending(bool flag) { fPending = flag; }

  /**
   * @brief Function that set next path
   *
   * @param path Volume path
   */
  VECCORE_ATT_HOST_DEVICE
  void SetPath(VolumePath_t const *const path);

  /**
   * @brief Function that set next volume path
   *
   * @param path Volume path
   */
  VECCORE_ATT_HOST_DEVICE
  void SetNextPath(VolumePath_t const *const path);

  /** @brief return the contiguous memory size needed to hold a GeantTrack*/
  VECCORE_ATT_HOST_DEVICE
  static size_t SizeOfInstance();

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack *MakeInstanceAt(void *addr);

  /**
   * @brief GeantTrack MakeInstance allocating the necessary buffer on the heap.
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack *MakeInstance();

  /**
   * @brief Releases track instance created using MakeInstance.
   */
  VECCORE_ATT_HOST_DEVICE
  static void ReleaseInstance(GeantTrack *track);

  /**
   * @brief Rounds up an address to the aligned value
   * @param buf Buffer address to align
   */
  VECCORE_ATT_HOST_DEVICE
  static char *round_up_align(char *buf) {
    long remainder = ((long)buf) % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return buf;
    return (buf + GEANT_ALIGN_PADDING - remainder);
  }

  /**
   * @brief Rounds up a value to upper aligned version
   * @param buf Buffer address to align
   */
  VECCORE_ATT_HOST_DEVICE
  static size_t round_up_align(size_t value) {
    size_t remainder = ((size_t)value) % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return value;
    return (value + GEANT_ALIGN_PADDING - remainder);
  }
  //  ClassDefNV(GeantTrack, 1) // The track
};

// Track token for user data
class TrackToken
{
private:
  std::string fName;      /** token name */
  size_t fOffset;         /** offset with respect to track user data base address */
public:
  TrackToken(const char *name = "", size_t offset = 0) : fName(name), fOffset(offset) {}

  std::string const GetName() const { return fName; }

  template <typename T>
  T &Data(GeantTrack *track) {
    return *(T*)(track->fExtraData + fOffset);
  }
};

// A structure keeping a lambda function to construct user data in place, at
// a given offset with respect to a base address
struct DataInplaceConstructor_t {
  size_t fDataOffset = 0;
  InplaceConstructFunc_t fDataConstructor;
  PrintUserDataFunc_t fPrintUserData;

  void operator()(GeantTrack &track) {
    fDataConstructor(track.fExtraData + fDataOffset);
  }

  void PrintUserData(GeantTrack const &track) const {
    fPrintUserData(track.fExtraData + fDataOffset);
  }

  DataInplaceConstructor_t(size_t offset, InplaceConstructFunc_t ctor, PrintUserDataFunc_t prfunc)
    : fDataOffset(offset), fDataConstructor(ctor), fPrintUserData(prfunc) {}
};

// User data manager singleton. Becomes read-only after registration phase
class TrackDataMgr {
private:
  bool fLocked = false;   /** Lock flag. User not allowed to register new data type.*/
  size_t fTrackSize = 0;  /** Total track size including alignment rounding */
  size_t fDataOffset = 0; /** Offset for the user data start address with respect to GeantTrack::fExtraData pointer */
  size_t fMaxDepth = 0;   /** Maximum geometry depth, to be provided at construction time */

  vector_t<DataInplaceConstructor_t> fConstructors; /** Inplace user-defined constructors to be called at track initialization. */
  vector_t<TrackToken> fTokens;
  std::mutex fRegisterLock; /** Multithreading lock for the registration phase */

  VECCORE_ATT_HOST_DEVICE
  TrackDataMgr(size_t maxdepth);

public:
  VECCORE_ATT_HOST_DEVICE
  static
  TrackDataMgr *GetInstance(size_t fMaxDepth = 0);

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  size_t GetMaxDepth() const { return fMaxDepth; }

  GEANT_FORCE_INLINE
  size_t GetTrackSize() const { return fTrackSize; }

  GEANT_FORCE_INLINE
  size_t GetDataSize() const { return fDataOffset; }

  GEANT_FORCE_INLINE
  const TrackToken *GetToken(const char *name) const
  {
    for (size_t i=0; i<fTokens.size(); ++i) { if (fTokens[i].GetName() == name) return &fTokens[i]; }
    return nullptr;
  }

  /** @brief Apply in-place construction of user-defined data on the track */
  void InitializeTrack(GeantTrack &track) {
    fLocked = true;
    // Construct objects at base address + offset
    for (auto construct : fConstructors)
      construct(track);
  }

  /** @brief Invoke all registered print methods for user data */
  void PrintUserData(GeantTrack const &track) const {
    for (auto userdata : fConstructors)
      userdata.PrintUserData(track);
  }

  /** @brief Lock any further registration */
  void Lock()
  {
    std::lock_guard<std::mutex> lock(fRegisterLock);
    fLocked = true;
  }

  void Print()
  {
    std::cout << "*** TrackDataMgr report: track size = " << fTrackSize << " bytes,  max. depth = " << fMaxDepth << std::endl;
    std::cout << "                         extra data size = " <<  fDataOffset << " bytes in the following blocks: ";
    for (const auto &token : fTokens) std::cout << token.GetName() << " ";
    std::cout << "\n";
  }

  template <typename T /*, typename... Params*/>
  TrackToken RegisterDataType(const char *name /*, Params... params*/)
  {
    if (fLocked) throw std::runtime_error("Registering types in the track is locked");
    std::lock_guard<std::mutex> lock(fRegisterLock);

    // If a token with same name is found, return its reference
    const TrackToken *token = GetToken(name);
    if (token) return (*token);

    // Register lambda for constructing in-place the user data
    auto userinplaceconstructor = [&](void * address) {
      new (address) T(/*params...*/);
    };

    // Register print function for debugging
    auto printuserdata = [&](void * address) {
      T const &userdata = *(T*)(address);
      userdata.Print();
    };

    // Calculate data offset
    size_t block_size = GeantTrack::round_up_align(sizeof(T)); // multiple of the alignment padding
    fTrackSize += block_size;
    DataInplaceConstructor_t ctor {fDataOffset, userinplaceconstructor, printuserdata};
    fConstructors.push_back(ctor);

    // Create a new token
    TrackToken newtoken(name, fDataOffset);
    fTokens.push_back(newtoken);
    fDataOffset += block_size;
    return (newtoken);
  }
};

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
