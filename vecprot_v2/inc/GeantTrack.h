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

#include "Geant/Typedefs.h"

#include <atomic>

#ifndef GEANT_ALIGN_PADDING
#define GEANT_ALIGN_PADDING 32
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

const double kB2C = -0.299792458e-3;

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

class GeantTrack_v;
class GeantTaskData;

/*
#if ROOT_VERSION_CODE <= ROOT_VERSION(6, 04, 00)
// Work around bug in ROOT v6.04
#ifdef G__DICTIONARY
class GeantTrack;
using Track = GeantTrack;
using Track_v = GeantTrack_v;
#endif
#endif
*/

/**
 * @brief Class GeantTrack
 */
class GeantTrack {
public:
  int fEvent;            /** Event number */
  int fEvslot;           /** Event slot */
  int fParticle;         /** Index of corresponding particle */
  int fPDG;              /** Particle pdg code */
  int fGVcode;           /** GV particle code */
  int fEindex;           /** Element index */
  int fCharge;           /** Particle charge */
  int fProcess;          /** Current process */
  int fVindex;           /** Current volume index */
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
  double Beta() const { return fP / fE; }

  /** @brief Function that return charge value */
  int Charge() const { return fCharge; }

  /** @brief Function that return curvature */
  double Curvature() const;

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
  double Gamma() const { return fMass ? fE / fMass : std::numeric_limits<double>::max(); }

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

  /** @brief Function that return current volume index */
  int Vindex() const { return fVindex; }

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
   * @brief Function that read track from vector
   *
   * @param arr Array of tracks
   * @param i Position to read
   */
  void ReadFromVector(const GeantTrack_v &arr, int i);

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
   * @brief Function that set current volume index
   *
   * @param ind Current volume index that should be set as fVindex
   */
  void SetVindex(int ind) { fVindex = ind; }

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

/**
 * @brief SOA for GeantTrack used at processing time
 *
 */
class GeantTrack_v {
public:
  static size_t const cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];
#ifdef GEANT_CUDA_DEVICE_BUILD
  int fNtracks; /** Number of tracks contained */
#else
  std::atomic_int fNtracks; /** number of tracks contained */
#endif
  int fNselected; /** Number of selected tracks */
  bool fCompact;  /** Flag marking the compactness */
  bool fMixed;    /** Contains tracks in mixed volumes */
  cacheline_pad_t pad0_;
  int fMaxtracks;    /** Max size for tracks */
  BitSet *fHoles;    /** Bits of holes */
  BitSet *fSelected; /** Mask of selected tracks for the current operation */

  int fMaxDepth;   /** Maximum geometry depth allowed */
  size_t fBufSize; /** Size of the internal buffer */
  char *fVPstart;  /** Address of volume path buffer */
  char *fBuf;      /** Buffer holding tracks data */

  int *fEventV;            /** Event numbers */
  int *fEvslotV;           /** Event slots */
  int *fParticleV;         /** Indices of corresponding particles */
  int *fPDGV;              /** Particle pdg codes */
  int *fGVcodeV;           /** GV internal codes */
  int *fEindexV;           /** Element indices */
  int *fChargeV;           /** Particle charges */
  int *fProcessV;          /** Current process */
  int *fVindexV;           /** Volume index */
  int *fNstepsV;           /** Number of steps made */
  Species_t *fSpeciesV;    /** Particle species */
  TrackStatus_t *fStatusV; /** Track statuses */
  double *fMassV;          /** Particle masses */
  double *fXposV;          /** Array of track X positions */
  double *fYposV;          /** Array of track Y positions */
  double *fZposV;          /** Array of track Z positions */
  double *fXdirV;          /** Array of track directions */
  double *fYdirV;
  double *fZdirV;
  double *fPV;      /** Momentum */
  double *fEV;      /** Energies */
  double *fTimeV;   /** Time */
  double *fEdepV;   /** Energy depositions */
  double *fPstepV;  /** Selected physical steps */
  double *fStepV;   /** Current steps */
  double *fSnextV;  /** Straight distances to next boundary */
  double *fSafetyV; /** Safe distances to any boundary */
  double *fNintLenV;/** Number of interaction lenghts traveled in last step */
  double *fIntLenV; /** Cumulated interaction length since last discrete process */
  bool *fBoundaryV; /** True if starting from boundary */
  bool *fPendingV;
  VolumePath_t **fPathV;     /** Paths for the particles in the geometry */
  VolumePath_t **fNextpathV; /** Paths for next volumes */

  /**
   * @brief Function that assign in current buffer
   *
   * @param buff Buffer to be assigned to
   * @param size Size of buffer
   */
  GEANT_CUDA_BOTH_CODE
  void AssignInBuffer(char *buff, int size);

  /**
   * @brief Function that copy to current buffer
   *
   * @param buff Buffer to be copied to
   * @param size Size of buffer
   */
  void CopyToBuffer(char *buff, int size);

private:
  /** @brief Copy constructor (not allowed) */
  GeantTrack_v(const GeantTrack_v &track_v);

  /** @brief Operator = */
  GeantTrack_v &operator=(const GeantTrack_v &track_v);

  /**
   * @brief GeantTrack constructor based on a provided single buffer.
   */
  GEANT_CUDA_BOTH_CODE
  GeantTrack_v(void *addr, unsigned int nTracks, int maxdepth);

public:
  /** @brief GeantTrack_v constructor */
  GeantTrack_v();

  /**
   * @brief GeantTrack_v parametrized constructor
   *
   * @param size Size of track
   * @param maxdepth Maximum allowed geometry depth
   */
  GeantTrack_v(int size, int maxdepth);

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  GEANT_CUDA_BOTH_CODE
  static GeantTrack_v *MakeInstanceAt(void *addr, unsigned int nTracks, int maxdepth);

  /** @brief GeantTrack_v destructor */
  ~GeantTrack_v();

  /** @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth */
  GEANT_CUDA_BOTH_CODE
  static size_t SizeOfInstance(size_t nTracks, size_t maxdepth);

  /** @brief  Function that returns the buffer size  */
  size_t BufferSize() const { return fBufSize; }

  /** @brief  Return the address of the data memory buffer  */
  void *Buffer() const { return fBuf; }

  /** @brief  Function that returns buffer size needed to hold the data for nTracks and maxdepth */
  GEANT_CUDA_BOTH_CODE
  static size_t BufferSize(size_t nTracks, size_t maxdepth);

  /** @brief  Function that returned max size for tracks */
  int Capacity() const { return fMaxtracks; }

  /** @brief  Check direction normalization within tolerance */
  bool IsNormalized(int itr, double tolerance = 1.E-8) const;

  /**
   * @brief Function that compare 2 tracks
   *
   * @param tr1 First track
   * @param i1 Bit number 'i1'
   * @param tr2 Second track
   * @param i2 Bit number 'i2'
   */
  static bool IsSame(const GeantTrack_v &tr1, int i1, const GeantTrack_v &tr2, int i2);

#ifdef GEANT_CUDA_DEVICE_BUILD

  /** @brief  Function that returned number of tracks contained  */
  GEANT_CUDA_BOTH_CODE
  int GetNtracks() const { return fNtracks; }

  /** @brief  Function that set number of tracks contained  */
  GEANT_CUDA_BOTH_CODE
  void SetNtracks(int ntracks) { fNtracks = ntracks; }
#else

  /** @brief  Function that returned number of tracks contained  C++11 */
  int GetNtracks() const { return fNtracks.load(); }

  /** @brief  Function that set number of tracks contained  C++11 */
  void SetNtracks(int ntracks) { fNtracks.store(ntracks); }
#endif

  /** @brief  Function that return number of selected tracks  */
  int GetNselected() const { return fNselected; }

  /**
   * @brief Add track function
   *
   * @param track Track that should be added
   * @param import Flag for importing (by default False)
   */
  GEANT_CUDA_BOTH_CODE
  int AddTrack(GeantTrack &track, bool import = false);

  /**
   * @brief Add & synctrack function
   *
   * @param track Track that should be added
   */
  int AddTrackSync(GeantTrack &track);

  /**
   * @brief Add track function
   *
   * @param arr Array of tracks
   * @param i  Bit number 'i'
   * @param import Flag for importing (by default False)
   */
  GEANT_CUDA_BOTH_CODE
  int AddTrack(GeantTrack_v &arr, int i, bool import = false);

  /**
   * @brief Add & synctrack function
   *
   * @param arr Track array
   * @param i Bit number 'i'
   */
  int AddTrackSync(GeantTrack_v &arr, int i);

  /**
   * @brief AddAt & synctrack function
   *
   * @param itrack which location to write the track to
   * @param arr input Track array
   * @param i Bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  int AddTrackSyncAt(int itrack, GeantTrack_v &arr, int i);

  /**
   * @brief Add track function
   *
   * @param arr Track array
   * @param istart Start index of tracks (Start bit number 'i')
   * @param iend End index of tracks (End bit number 'i')
   * @param import Import flag (by default False)
   */
  void AddTracks(GeantTrack_v &arr, int istart, int iend, bool import = false);

  /** @brief Function that check track */
  void CheckTracks();

  /**
   * @brief  Function that mark removed bit number 'i' through value of bits of holes and flag marking the compactness
   *
   * @param i Bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  void MarkRemoved(int i) {
    fHoles->SetBitNumber(i);
    fCompact = false;
  }

  /**
   * @brief Function that removes tracks
   *
   * @param from Start index of tracks to be removed
   * @param to End index of tracks to be removed
   */
  void RemoveTracks(int from, int to);

  /** @brief Function that return size of track */
  size_t Sizeof() const { return sizeof(GeantTrack_v) + fBufSize; }

  /**
   * @brief Function that delete track
   *
   * @param itr Track ID
   */
  void DeleteTrack(int itr);

  /**
   * @brief Function that deselect bit number 'i' from value
   *
   * @param i Bit number i
   */
  void Deselect(int i) { fSelected->SetBitNumber(i, false); }

  /** @brief Function that deselect all bit number from value  */
  void DeselectAll() {
    fSelected->ResetAllBits();
    fNselected = 0;
  }

  /** @brief Function to normalize direction */
  void Normalize(int itr) __attribute__((always_inline)) {
    double norm = 1. / Math::Sqrt(fXdirV[itr] * fXdirV[itr] + fYdirV[itr] * fYdirV[itr] + fZdirV[itr] * fZdirV[itr]);
    fXdirV[itr] *= norm;
    fYdirV[itr] *= norm;
    fZdirV[itr] *= norm;
  }

  /**
   * @brief Function that select bit number 'i' to be value
   *
   * @param i Bit number i
   */
  void Select(int i) { fSelected->SetBitNumber(i); }

  /**
   * @brief  Function that select tracks
   * @details Function that select tracks
   *
   * @param n Number of tracks to be selected
   */
  void SelectTracks(int n) { fNselected = n; }

  /**
   * @brief Set if tracks are in mixed volumes through fMixed flag
   *
   * @param flag Flag that shows if it is contained tracks in mixed volumes
   */
  void SetMixed(bool flag) { fMixed = flag; }

  /**
   * @brief Sorting function for track status
   *
   * @param status Track status to be checked
   */
  int SortByStatus(TrackStatus_t status);

  /**
   * @brief Sorting function for tracks where the step was limited by discrete processes
   *
   * @return Number of selected tracks
   */
  int SortByLimitingDiscreteProcess();
  /**
   * @brief Function for removal tracks according status
   *
   * @param status Track status to be selected for removal
   * @param output Output array of tracks
   */
  int RemoveByStatus(TrackStatus_t status, GeantTrack_v &output);

  /**
   * @brief Selection function
   *
   * @param i Bit number 'i' to be value
   */
  bool IsSelected(int i) { return fSelected->TestBitNumber(i); }

  /** @brief Clear function */
  GEANT_CUDA_BOTH_CODE
  void Clear(const char *option = "");

  /**
   * @brief Compact function
   * @param  moveto ???????
   */
  int Compact(GeantTrack_v *moveto = 0);

  /**
   * @brief Contain function
   *
   * @param evstart Event to start from
   * @param nevents Numbe rof evets (by default 1)
   */
  bool Contains(int evstart, int nevents = 1) const;

  /** @brief Clear selection of bit value */
  void ClearSelection() { fSelected->ResetAllBits(); }

  /**
   * @brief Function that return track
   *
   * @param i Bit number 'i'
   * @param track Track to be returned
   */
  void GetTrack(int i, GeantTrack &track) const;

  /** @brief Function that check flag marking the compactness */
  bool IsCompact() const { return fCompact; }

  /** @brief Function that return how many tracks in mixed volumes */
  bool IsMixed() const { return fMixed; }

  /** @brief Function that print pointers */
  void PrintPointers() { printf("fEventV=%p fBoundaryV=%p\n", (void *)fEventV, (void *)fBoundaryV); }

  /**
   * @brief Function that print track
   *
   * @param itr Track ID
   */
  GEANT_CUDA_BOTH_CODE
  void PrintTrack(int itr, const char *msg = "") const;

  /** @brief Function that print all tracks */
  void PrintTracks(const char *msg = "") const;

  /**
   * @brief Function for navigation that check if location is the same or not for single track
   *
   * @param itr Track ID
   * @param start Start volume path
   * @param end End volume path
   */

// void InspectGeometryState(int itr) const;
// void InspectIsSameLocation(int itr) const;

#ifdef USE_VECGEOM_NAVIGATOR

  /**
   * @brief Function that check if location path's consistency
   *
   * @param itr Track ID
   */
  void CheckLocationPathConsistency(int itr) const;
#endif

  /**
   * @brief Function that provides postponed action for tracks
   *
   * @param ntracks Number of tracks
   */
  TransportAction_t PostponedAction(int ntracks) const;

  /**
   * @brief Function that provides postponed action for track
   *
   * @param itr Track ID
   * @param output Output of tracks
   */
  GEANT_CUDA_BOTH_CODE
  int PostponeTrack(int itr, GeantTrack_v &output);

  /**
   * @brief Function that provides postponed action for all track
   *
   * @param output Output of tracks
   */
  int PostponeTracks(GeantTrack_v &output);
  // void      PropagateBack(int itr, double crtstep);

  /**
   * @brief Function that compute transport length
   *
   * @param ntracks Number of tracks and TaskData object ( with preallocated thread/task local workspaces )
   */
  GEANT_CUDA_BOTH_CODE
  void ComputeTransportLength(int ntracks, GeantTaskData *);

  /**
   * @brief Function that compute single transport length ?????
   *
   * @param itr Track ID
   */
  GEANT_CUDA_BOTH_CODE
  void ComputeTransportLengthSingle(int itr, GeantTaskData *);

  /**
   * @brief Function of propagation in volume
   *
   * @param ntracks Number of tracks
   * @param crtstep ??????
   * @param tid Track ID
   */
  GEANT_CUDA_BOTH_CODE
  void PropagateInVolume(int ntracks, const double *crtstep, GeantTaskData *td);

  /**
   * @brief Function of propagation of track in volume
   *
   * @param i Bit number 'i'
   * @param crtstep ???????
   * @param tid Track ID
   */
  GEANT_CUDA_BOTH_CODE
  void PropagateInVolumeSingle(int i, double crtstep, GeantTaskData *td);

  /**
   * @brief Popagation function in straight trajectories
   *
   * @param ntracks Number of tracks
   * @param crtstep  ?????
   */
  int PropagateStraight(int ntracks, double *crtstep);

  /**
   * @brief Function of propagation of tracks
   *
   * @param output Output array of tracks
   * @param tid Track ID
   */
  int PropagateTracks(GeantTaskData *td);

  GEANT_CUDA_BOTH_CODE
  int PropagateTracksScalar(GeantTaskData *td, int stage = 0);

  GEANT_CUDA_BOTH_CODE
  int PropagateSingleTrack(int itr, GeantTaskData *td, int stage);

  /**
   * @brief Resize function
   * @param newsize New size to be resized
   */
  void Resize(int newsize);

  /** @brief Function that replace track at positions i with track on position j */
  void ReplaceTrack(int i, int withj);

  /** @brief Reshuffle function */
  int Reshuffle();

  /**
   * @brief Function that swap tracks at positions i and j
   * @param i Input bit number 'i'
   * @param j Input bit number 'j'
   */
  void SwapTracks(int i, int j);

  /**
   * @brief Function that return beta value
   * @param  i Input bit number 'i'
   */
  double Beta(int i) const { return fPV[i] / fEV[i]; }

  /**
   * @brief Function that return curvature in different areas of geometry
   * @param  i Input bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  double Curvature(int i) const;

  /** @brief Function that return safe length */
  GEANT_CUDA_BOTH_CODE
  double SafeLength(int i, double eps = 1.E-4);

  /**
   * @brief Function that return gamma value
   * @param  i Input bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  double Gamma(int i) const { return fMassV[i] ? fEV[i] / fMassV[i] : std::numeric_limits<double>::max(); }

  /**
   * @brief Function that return X projection of momentum value
   * @param  i Input bit number 'i'
   */
  double Px(int i) const { return fPV[i] * fXdirV[i]; }

  /**
   * @brief Function that return Y projection of momentum value
   * @param  i Input bit number 'i'
   */
  double Py(int i) const { return fPV[i] * fYdirV[i]; }

  /**
   * @brief Function that return Z projection of momentum value
   * @param  i Input bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  double Pz(int i) const { return fPV[i] * fZdirV[i]; }

  /**
   * @brief Function that return module of momentum value
   * @param  i Input bit number 'i'
   */
  GEANT_CUDA_BOTH_CODE
  double Pt(int i) const { return fPV[i] * Math::Sqrt(fXdirV[i] * fXdirV[i] + fYdirV[i] * fYdirV[i]); }

  /** @brief Function that return time traveled in the step */
  GEANT_CUDA_BOTH_CODE
  double TimeStep(int i, double step) const { return fEV[i]*step/fPV[i]; }

  /**
   * @brief Function that returns the logical volume of the i-th track
   * @param  i Input bit number 'i'
   */
  Volume_t const*GetVolume(int i) const;

  /**
   * @brief Function that returns next logical volume of i-th track
   * @param  i Input bit number 'i'
   */
  Volume_t const*GetNextVolume(int i) const;

  /**
   * @brief Function that returns the current material the i-th track is in
   * @param  i Input bit number 'i'
   */
  Material_t *GetMaterial(int i) const;

  /** @brief Function allowing to set a breakpoint on a given step */
  GEANT_CUDA_BOTH_CODE
  bool BreakOnStep(int evt, int trk, int stp, int nsteps = 1, const char *msg = "", int itr = -1);

  /**
   * @brief Check consistency of track navigation
   * @param  itr Track number to be checked
   */
  bool CheckNavConsistency(int itr);

  /**
   * @brief Function round up align ?????
   * @param num Number ?????
   */
  GEANT_CUDA_BOTH_CODE
  static int round_up_align(int num) {
    int remainder = num % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return num;
    return (num + GEANT_ALIGN_PADDING - remainder);
  }

  GEANT_CUDA_BOTH_CODE
  static char *round_up_align(char *buf) {
    long remainder = ((long)buf) % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return buf;
    return (buf + GEANT_ALIGN_PADDING - remainder);
  }

  //  ClassDefNV(GeantTrack_v, 1) // SOA for GeantTrack class
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA
#ifdef GEANT_NVCC
namespace cxx {
class GeantTrack_v;
}
#else
namespace cuda {
class GeantTrack_v;
}
#endif

bool ToDevice(vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> dest, cxx::GeantTrack_v *source, cudaStream_t stream);
bool FromDevice(cxx::GeantTrack_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> source, cudaStream_t stream);
void FromDeviceConversion(cxx::GeantTrack_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrack_v> source);
#endif

} // Geant

#endif
