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
#include <VecCore/VecCore>
#include "Geant/Config.h"
#include "Geant/math_wrappers.h"

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
  kAtRestActionStage,       // At-rest action stage (at-rest part of the inetraction)
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

constexpr size_t kNstages = size_t(kSteppingActionsStage) + 1;
constexpr size_t kNumPhysicsProcess = 10;

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

template <typename T>
using Vector3D = vecgeom::Vector3D<T>;

private:
  int fEvent = -1;           /** Event number */
  int fEvslot = -1;          /** Event slot */
  int fParticle = -1;        /** Index of corresponding particle */
  int fPrimaryIndx = -1;     /** Index of the primary particle in the current event */
  int fMother = -1;           /** Index of mother particle */
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
  bool fIsOnBoundaryPreStp = false;  /** to indicate that the particle was on boundary at the pre-step pint */
  bool fPrePropagationDone = false;  /** Indicate if pre-propagation stage was done for this particle. */
  Volume_t const *fVolume = nullptr; /** Current volume the particle is in */

  // max number of physics processesper particle is assumed to be 10!
  int  fPhysicsProcessIndex = -1;  // selected physics process
  double  fPhysicsNumOfInteractLengthLeft[kNumPhysicsProcess];
  double  fPhysicsInteractLength[kNumPhysicsProcess]; // mfp

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

  static constexpr double kB2C = -0.299792458e-3;
  static constexpr double kTiny = 1.E-50;

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

  /** @brief GeantTrack copy constructor */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(const GeantTrack &other) = delete;

  /** @brief Operator = */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack &operator=(const GeantTrack &other);

  //---------------------------------//
  //*** Basic and derived getters ***//
  //---------------------------------//

  /** @brief Getter for the event number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Event() const { return fEvent; }

  /** @brief Getter for the slot number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int EventSlot() const { return fEvslot; }

  /** @brief Getter for the index of corresponding particle */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Particle() const { return fParticle; }

  /** @brief Getter for the index of the primary particle in the current event */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int PrimaryParticleIndex() const { return fPrimaryIndx; }

  /** @brief Getter for thes index of mother particle */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Mother() const { return fMother; }

  /** @brief Getter for the particle pdg code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int PDG() const { return fPDG; }

  /** @brief Getter for the GV particle code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GVcode() const { return fGVcode; }

  /** @brief Getter for the element index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int EIndex() const { return fEindex; }

  /** @brief Getter for the index in the track block */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int BIndex() const { return fBindex; }

  /** @brief Getter for the charge value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Charge() const { return fCharge; }

  /** @brief Fills scalar charge from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD charges to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetCharge_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Real_v &charge_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(charge_v, i, tracks[offset + i]->Charge());
    return nelem;
  }

  /** @brief Getter for the current process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int Process() const { return fProcess; }

  /** @brief Getter for the number of physical step made */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetNsteps() const { return fNsteps; }

  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetMaxDepth() const { return fMaxDepth; }

  /** @brief Getter for simulation stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  ESimulationStage Stage() const { return (ESimulationStage)fStage; }

  /** @brief Getter for simulation stage as integer */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetStage() const { return fStage; }

  /** @brief Getter for track generation (0 = primary)*/
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetGeneration() const { return fGeneration; }

  /** Getter for the particle species */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Species_t Species() const { return fSpecies; }

  /** Getter for the track status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  TrackStatus_t Status() const { return fStatus; }

  /** @brief Getter for the rest mass value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Mass() const { return fMass; }

  /** @brief X position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double X() const { return fXpos; }

  /** @brief Y position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Y() const { return fYpos; }

  /** @brief Z position */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Z() const { return fZpos; }

  /** @brief Getter for the pointer to X position value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  const double *Position() const { return &fXpos; }

  /** @brief Fills scalar position components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param pos_v SIMD position to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetPos_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Vector3D<Real_v> &pos_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(pos_v[0], i, tracks[offset + i]->X());
      vecCore::Set(pos_v[1], i, tracks[offset + i]->Y());
      vecCore::Set(pos_v[2], i, tracks[offset + i]->Z());
    }
    return nelem;
  }

  /** @brief Getter for the X direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Dx() const { return fXdir; }

  /** @brief Getter for the Y direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Dy() const { return fYdir; }

  /** @brief Getter for the Z direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Dz() const { return fZdir; }

  /** @brief Getter for the pointer to X direction value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  const double *Direction() const { return &fXdir; }

  /** @brief Fills scalar direction components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param dir_v SIMD direction to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetDir_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Vector3D<Real_v> &dir_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(dir_v[0], i, tracks[offset + i]->Dx());
      vecCore::Set(dir_v[1], i, tracks[offset + i]->Dy());
      vecCore::Set(dir_v[2], i, tracks[offset + i]->Dz());
    }
    return nelem;
  }

  /** @brief Getter for the momentum value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double P() const { return fP; }

  /** @brief Getter for the momentum X component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Px() const { return fP * fXdir; }

  /** @brief Getter for the momentum Y component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Py() const { return fP * fYdir; }

  /** @brief Getter for the momentum Z component */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Pz() const { return fP * fZdir; }

  /** @brief Fills scalar momentum components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param pos_v SIMD momentum to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetP_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Vector3D<Real_v> &mom_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i) {
      vecCore::Set(mom_v[0], i, tracks[offset + i]->Px());
      vecCore::Set(mom_v[1], i, tracks[offset + i]->Py());
      vecCore::Set(mom_v[2], i, tracks[offset + i]->Pz());
    }
    return nelem;
  }

  /** @brief Getter for the module momentum's value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Pt() const { return fP * Math::Sqrt(fXdir * fXdir + fYdir * fYdir); }

  /** @brief Getter for the energy value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double E() const { return fE; }

  /** @brief Fills scalar energy components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD energies to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetE_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Real_v &e_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(e_v, i, tracks[offset + i]->E());
    return nelem;
  }

  /** @brief Getter for the kinetic energy value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double T() const { return (fE - fMass); }

  /** @brief Fills scalar kinetic energy components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param t_v SIMD kinetic energies to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetT_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Real_v &t_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(t_v, i, tracks[offset + i]->T());
    return nelem;
  }

  /** Getter for the time */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Time() const { return fTime; }

  /** @brief Getter for the energy deposition value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Edep() const { return fEdep; }

  /** @brief Getter for the selected physical step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetPstep() const { return fPstep; }

  /** @brief Getter for the physical step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetStep() const { return fStep; }

  /** @brief Getter for the time traveled in the step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double TimeStep(double step) const { return fE*step/fP; }

  /** @brief Getter for the straight distance to next boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetSnext() const { return fSnext; }

  /** @brief Getter for the safe distance to any boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetSafety() const { return fSafety; }

  /** @brief Getter for the number of interaction lengths travelled since last step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetNintLen() const { return fNintLen; }

  /** @brief Getter for the interaction length since last discrete process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetIntLen() const { return fIntLen; }

  /** @brief Getter for the true if starting from boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool Boundary() const { return fBoundary; }

  /** @brief Getter for the pending track status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool Pending() const { return fPending; }

  /** @brief Getter for path ownership */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsOwnPath() const { return fOwnPath; }

  /** @brief Getter pre step boundary status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsOnBoundaryPreStp() const { return fIsOnBoundaryPreStp; }

  /** @brief Getter for pre propagation done status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  bool IsPrePropagationDone() const { return fPrePropagationDone; }

  /** @brief Getter for the volume */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  Volume_t const *GetVolume() const { return fVolume; }

  /** @brief Getter for the physics process index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int GetPhysicsProcessIndex() const { return fPhysicsProcessIndex; }

  /** @brief Get number of remaining interaction lengths to the physics process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetPhysicsNumOfInteractLengthLeft(int iproc) const { return fPhysicsNumOfInteractLengthLeft[iproc]; }

  /** @brief Get total interaction length to the physics process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double GetPhysicsInteractLength(int iproc) const { return fPhysicsInteractLength[iproc]; }

  /** @brief Getter for thes current path */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  VolumePath_t *Path() const { return fPath; }

  /** @brief Getter for thes next path */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  VolumePath_t *NextPath() const { return fNextpath; }

  /** @brief Getter for the next volume (no backup) */
  VECCORE_ATT_HOST_DEVICE
  Volume_t const *GetNextVolume() const {
#ifdef USE_VECGEOM_NAVIGATOR
    return ( fNextpath->Top()->GetLogicalVolume() );
#else
    return ( fNextpath->GetCurrentNode()->GetVolume() );
#endif
  }

  /** @brief Getter for the material */
  VECCORE_ATT_HOST_DEVICE
  Material_t *GetMaterial() const;

  /** @brief Getter for the beta value */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Beta() const { return fP / fE; }

  /** @brief Fills scalar beta from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param beta_v SIMD betas to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetBeta_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Real_v &beta_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(beta_v, i, tracks[offset + i]->Beta());
    return nelem;
  }

  /** @brief Getter for the curvature. To be changed when handling properly field*/
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Curvature(double Bz) const {
    // Curvature
    double qB = fCharge * Bz;
    if (fabs(qB) < kTiny) return kTiny;
    return fabs(kB2C * qB / (Pt() + kTiny));
  }

  /** @brief Getter for the gamma value*/
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  double Gamma() const { return fMass ? fE / fMass : DBL_MAX; }

  /** @brief Fills scalar gamma components from input vector into SIMD type.
   *  @param tracks Vector of pointers to tracks
   *  @param offset Start offset.
   *  @param ntracks Number of tracks in the input vector
   *  @param e_v SIMD gammas to fill
   *  @return Number of filled lanes. User should check if deciding it is worth to extract the tail.
   **/
  template <typename Real_v, bool Tail = false>
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  static size_t GetGamma_v(TrackVec_t const &tracks, const size_t offset, const size_t ntracks, Real_v &gamma_v)
  {
    constexpr size_t kVecSize = vecCore::VectorSize<Real_v>();
    size_t nelem = kVecSize;
    if (Tail) nelem = ntracks - offset;
    assert(offset <= ntracks - nelem );
    for (size_t i = 0; i < nelem; ++i)
      vecCore::Set(gamma_v, i, tracks[offset + i]->Gamma());
    return nelem;
  }

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

  //---------------------//
  //*** Basic setters ***//
  //---------------------//

  /** @brief Setter for event number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEvent(int event) { fEvent = event; }

  /** @brief Setter for event slot number */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEvslot(int slot) { fEvslot = slot; }

  /** @brief Setter for particle id */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetParticle(int particle) { fParticle = particle; }

  /** @brief Setter for primary particle index in the current event */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPrimaryParticleIndex(int primaryindx) { fPrimaryIndx = primaryindx; }

  /** @brief Setter for mother index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetMother(int mother) { fMother = mother; }

  /** @brief Setter for particle pdg code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPDG(int pdg) { fPDG = pdg; }

  /** @brief Setter for GV particle code */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetGVcode(int gVcode) { fGVcode = gVcode; }

  /** @brief Setter for element index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEindex(int ind) { fEindex = ind; }

  /** @brief Setter for index of the track block */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBindex(int ind) { fBindex = ind; }

  /** @brief Setter for charge */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetCharge(int charge) { fCharge = charge; }

  /** @brief Setter for process index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetProcess(int process) { fProcess = process; }

  /** @brief Setter for number of steps */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNsteps(int nsteps) { fNsteps = nsteps; }

  /** @brief Increment number of steps */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void IncrementNsteps(int nsteps = 1) { fNsteps += nsteps; }

  /** @brief Setter for stage */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStage(int stage) { fStage = stage; }

  /** @brief Setter for particle generation */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetGeneration(int generation) { fGeneration = generation; }

  /** @brief Setter for particle species */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSpecies(Species_t species) { fSpecies = species; }

  /** @brief Setter for track status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStatus(TrackStatus_t status) { fStatus = status; }

  /** @brief Setter for particle mass */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetMass(double mass) { fMass = mass; }

  /** @brief Setter for position from components */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPosition(double x, double y, double z) {
    fXpos = x;
    fYpos = y;
    fZpos = z;
  }

  /** @brief Setter for position from vector */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPosition(Vector3D<double> const &pos) {
    fXpos = pos.x();
    fYpos = pos.y();
    fZpos = pos.z();
  }

  /** @brief Setter for direction from components */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetDirection(double dx, double dy, double dz) {
    fXdir = dx;
    fYdir = dy;
    fZdir = dz;
  }

  /** @brief Setter for direction from components */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetDirection(Vector3D<double> const &dir) {
    fXdir = dir.x();
    fYdir = dir.y();
    fZdir = dir.z();
  }

  /** @brief Setter for momentum */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetP(double p) { fP = p; }

  /** @brief Setter for energy */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetE(double e) { fE = e; }

  /** @brief Setter for energy */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void DecreaseE(double e) { fE -= e; }

  /** @brief Setter for time */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetTime(double time) { fTime = time; }

  /** @brief Increase time */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void IncreaseTime(double time) { fTime += time; }

  /** @brief Setter for energy deposition */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetEdep(double edep) { fEdep = edep; }

  /** @brief Setter for energy deposition */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void IncreaseEdep(double edep) { fEdep += edep; }

  /** @brief Setter for proposed physics step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPstep(double pstep) { fPstep = pstep; }

  /** @brief Decrease proposed physics step after some propagation */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void DecreasePstep(double len) { fPstep -= len; }

  /** @brief Setter for the physics step */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetStep(double step) { fStep = step; }

  /** @brief Increase progression step after some propagation */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void IncreaseStep(double len) { fStep += len; }

  /** @brief Setter for straight distance to next boundary */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSnext(double snext) { fSnext = snext; }

  /** @brief Decrease snext after some propagation */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void DecreaseSnext(double len) { fSnext -= len; }

  /** @brief Setter for the isotropic safe distance to volumes */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetSafety(double safety) { fSafety = safety; }

  /** @brief Decrease safety after some propagation */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void DecreaseSafety(double len) { fSafety -= len; }

  /** @brief Setter for the number of interaction lengths traveled */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNintLen(double nradlen) { fNintLen = nradlen; }

  /** @brief Setter for the interaction length since last discrete process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetIntLen(double intlen) { fIntLen = intlen; }

  /** @brief Setter for the starting from boundary flag */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBoundary(bool flag) { fBoundary = flag; }

  /** @brief Setter for the pending status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPending(bool flag) { fPending = flag; }

  /** @brief Setter for the on boundary pre-step status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetBoundaryPreStep(bool flag) { fIsOnBoundaryPreStp = flag; }

  /** @brief Setter for the pre propagation done status */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPrePropagationDone(bool flag) { fPrePropagationDone = flag; }

  /** @brief Setter for the physics process index */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPhysicsProcessIndex(int index) { fPhysicsProcessIndex = index; }

  /** @brief Setter for the number of physics interaction lengths left for a process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPhysicsNumOfInteractLengthLeft(int iproc, double val) { fPhysicsNumOfInteractLengthLeft[iproc] = val; }

  /** @brief Setter for the number of physics interaction lengths left for a process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void DecreasePhysicsNumOfInteractLengthLeft(int iproc, double val) { fPhysicsNumOfInteractLengthLeft[iproc] -= val; }

  /** @brief Setter for mean free path for a process */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetPhysicsInteractLength(int iproc, double val) { fPhysicsInteractLength[iproc] = val; }

  /** @brief Setter for the current geometry path */
  VECCORE_ATT_HOST_DEVICE
  void SetPath(VolumePath_t const *const path);

  /** @brief Setter for the next geometry path */
  VECCORE_ATT_HOST_DEVICE
  void SetNextPath(VolumePath_t const *const path);

  //---------------------//
  //***    Actions    ***//
  //---------------------//

  /** @brief Clone this track using specific task data storage */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack *Clone(GeantTaskData *td);

  /** @brief Function that stops the track depositing its kinetic energy */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Stop() { fEdep = T(); fE = fMass; fP = 0; }

  /** @brief Setter for the status killed to track */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Kill() { fStatus = kKilled; }

   /** Clear function */
  VECCORE_ATT_HOST_DEVICE
  void Clear(const char *option = "");

  /** Fast reset function */
  VECCORE_ATT_HOST_DEVICE
  void Reset(GeantTrack const &blueprint);

  /** @brief Print function */
  void Print(const char *msg = "") const;

  /** @brief Print function for a container of tracks */
  static void PrintTracks(TrackVec_t &tracks);

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

  /** @brief Returns the contiguous memory size needed to hold a GeantTrack*/
  VECCORE_ATT_HOST_DEVICE
  static size_t SizeOfInstance();

  /** @brief GeantTrack MakeInstance based on a provided single buffer */
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack *MakeInstanceAt(void *addr);

  /** @brief GeantTrack MakeInstance allocating the necessary buffer on the heap */
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack *MakeInstance();

  /** @brief Releases track instance created using MakeInstance */
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
