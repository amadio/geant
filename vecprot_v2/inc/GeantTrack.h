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
//  kGeometryStepStage,        // Compute geometry transport length
  kPrePropagationStage,      // Special msc stage for step limit phase
  kGeometryStepStage,        // Compute geometry transport length
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

GEANT_DECLARE_CONSTANT(double, gTolerance);

class GeantTaskData;
class GeantTrack;
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
  int fEvent;            /** Event number */
  int fEvslot;           /** Event slot */
  int fParticle;         /** Index of corresponding particle */
  int fMother;           /** Index of mother particle */
  int fPDG;              /** Particle pdg code */
  int fGVcode;           /** GV particle code */
  int fEindex;           /** Element index */
  int fBindex;           /** Index in the track block */
  int fCharge;           /** Particle charge */
  int fProcess;          /** Current process */
  int fNsteps;           /** Number of steps made */
  int fMaxDepth;         /** Maximum geometry depth */
  int fStage;            /** Simulation stage */
  int fGeneration;       /** Track generation: 0=primary */
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
  bool fPending;         /** Track pending to be processed */
  bool fOwnPath;         /** Marker for path ownership */
  VolumePath_t *fPath;   /** Paths for the particle in the geometry */
  VolumePath_t *fNextpath; /** Path for next volume */

// msc members: added here in mixed mode (int,bool,double) just for the development phase i.e. need to change later
//    -- NO SETTERS/GETTERS added
//----
//*****************
double fLambda0;        // elastic mean free path
double fLambda1;        // first transport mean free path
double fScrA;           // screening parameter if any
double fG1;             // first transport coef.
double fRange;
//******************
double fTheInitialRange;                // the initial (first step or first step in volume) range value of the particle
double fTheTrueStepLenght;              // the true step length
double fTheTransportDistance;           // the straight line distance between the pre- and true post-step points
double fTheZPathLenght;                 // projection of transport distance along the original direction
double fTheTrueGeomLimit;               // geometrical limit converted to true step length
double fTheDisplacementVectorX;         // displacement vector components X,Y,Z
double fTheDisplacementVectorY;
double fTheDisplacementVectorZ;
double fTheNewDirectionX;               // new direction components X,Y,Z (at the post-step point due to msc)
double fTheNewDirectionY;
double fTheNewDirectionZ;
double fPar1;
double fPar2;
double fPar3;

bool fIsEverythingWasDone;              // to indicate if everything could be done in the step limit phase
bool fIsMultipleSacettring;             // to indicate that msc needs to be perform (i.e. compute angular deflection)
bool fIsSingleScattering;               // to indicate that single scattering needs to be done
bool fIsEndedUpOnBoundary;              // ?? flag to indicate that geometry was the winer
bool fIsNoScatteringInMSC;              // to indicate that no scattering happend (i.e. forward) in msc
bool fIsNoDisplace;                     // to indicate that displacement is not computed
bool fIsInsideSkin;                     // to indicate that the particle is within skin from/to boundary
bool fIsWasOnBoundary;                  // to indicate that boundary crossing happend recently
bool fIsFirstStep;                      // to indicate that the first step is made with the particle
bool fIsFirstRealStep;                  // to indicate that the particle is making the first real step in the volume i.e.
                                        // just left the skin
//---


  /**
  * @brief GeantTrack in place constructor
  *
  * @param maxdepth Maximum geometry depth
  */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(void *addr, int maxdepth);

public:
  /** @brief GeantTrack constructor  */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack();

  /**
   * @brief GeantTrack copy constructor
   */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(const GeantTrack &other);

  /** @brief Operator = */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack &operator=(const GeantTrack &other);

  /**
  * @brief GeantTrack parametrized constructor
  *
  * @param ipdg PDG code
  * @param maxdepth Maximum geometry depth
  */
  VECCORE_ATT_HOST_DEVICE
  GeantTrack(int ipdg, int maxdepth);

  /** @brief GeantTrack destructor */
  VECCORE_ATT_HOST_DEVICE
  ~GeantTrack();

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
  Volume_t const*GetVolume() const;

  /** @brief Function that return next volume */
  VECCORE_ATT_HOST_DEVICE
  Volume_t const*GetNextVolume() const;

  /** @brief Function that return material */
  Material_t *GetMaterial() const;

  /** @brief Function that returns the current path (NavigationState) of the track */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  VolumePath_t *GetPath() const { return fPath; }

  /** @brief Function that return next path (NavigationState) of the track */
  VolumePath_t *GetNextPath() const { return fNextpath; }

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
   * @brief Function that set particle index
   *
   * @param particle Particle that should be set as fParticle
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetParticle(int particle) { fParticle = particle; }

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
   * @brief Function that set X, Y, Z component of the displacement vector provided by msc
   *
   * @param x X position
   * @param y Y position
   * @param z Z position
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetDisplacement(double x, double y, double z) {
    fTheDisplacementVectorX = x;
    fTheDisplacementVectorY = y;
    fTheDisplacementVectorZ = z;
  }


  /**
   * @brief Function that set the new X, Y, Z directions proposed by msc (will be applied or not depending other conditions)
   *
   * @param dx X direction
   * @param dy Y direction
   * @param dz Z direction
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void SetNewDirectionMsc(double dx, double dy, double dz) {
    fTheNewDirectionX = dx;
    fTheNewDirectionY = dy;
    fTheNewDirectionZ = dz;
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

  /** @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth */
  VECCORE_ATT_HOST_DEVICE
  static size_t SizeOfInstance(size_t maxdepth);

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack *MakeInstanceAt(void *addr, int maxdepth);

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

} // GEANT_IMPL_NAMESPACE
} // Geant

#endif
