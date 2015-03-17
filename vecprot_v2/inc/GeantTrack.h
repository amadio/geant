//===--- GeantTrack.h - Geant-V ---------------------------------*- C++ -*-===//
//
//                     Geant-V Prototype               
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTrack.h
 * @brief Implementation of track for Geant-V prototype 
 * @author Andrei Gheata 
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK
#define GEANT_TRACK

#include "Geant/Config.h"
#include "Geant/Math.h"

//#include "globals.h"
#include "TMath.h"
//#include "TBits.h"

#ifdef __STAT_DEBUG
#include "GeantTrackStat.h"
#endif

#if __cplusplus >= 201103L
#include <atomic>
#endif

#ifndef ALIGN_PADDING
#define ALIGN_PADDING 32
#endif

#ifndef VECCORE_BITSET_H
#include "BitSet.h"
typedef VecCore::BitSet BitSet;
#endif

#ifdef USE_VECGEOM_NAVIGATOR
#include "navigation/NavigationState.h"
typedef VECGEOM_NAMESPACE::NavigationState VolumePath_t;
#else
#include "TGeoBranchArray.h" // needed due to templated pools
typedef TGeoBranchArray VolumePath_t;
#endif

const Double_t kB2C = -0.299792458e-3;

/**
 * @enum TrackStatus_t
 */
enum TrackStatus_t {
  kAlive,
  kKilled,
  kInFlight,
  kBoundary,
  kExitingSetup,
  kPhysics,
  kPostponed,
  kNew
};

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

class TGeoMaterial;
class TGeoVolume;
class GeantTrack_v;

/**
 * @brief Class GeantTrack
 */
class GeantTrack {
public:
  Int_t fEvent;          /** Event number */
  Int_t fEvslot;         /** Event slot */
  Int_t fParticle;       /** Index of corresponding particle */
  Int_t fPDG;            /** Particle pdg code */
  Int_t fG5code;         /** G5 particle code */
  Int_t fEindex;         /** Element index */
  Int_t fCharge;         /** Particle charge */
  Int_t fProcess;        /** Current process */
  Int_t fVindex;         /** Current volume index */
  Int_t fNsteps;         /** Number of steps made */
  Species_t fSpecies;    /** Particle species */
  TrackStatus_t fStatus; /** Track status */
  Double_t fMass;        /** Particle mass */
  Double_t fXpos;        /** X position */
  Double_t fYpos;        /** Y position */
  Double_t fZpos;        /** Z position */
  Double_t fXdir;   /** X direction */
  Double_t fYdir;   /** Y direction */
  Double_t fZdir;   /** Z direction */
  Double_t fP;      /** Momentum */
  Double_t fE;      /** Energy */
  Double_t fTime;   /** Time */
  Double_t fEdep;   /** Energy deposition in the step */
  Double_t fPstep;  /** Selected physical step */
  Double_t fStep;   /** Current step */ 
  Double_t fSnext;  /** Straight distance to next boundary */
  Double_t fSafety; /** Safe distance to any boundary */
  Bool_t fFrombdr;  /** True if starting from boundary */
  Bool_t fPending;
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
  GeantTrack(Int_t ipdg);

  /** @brief GeantTrack destructor */
  ~GeantTrack();
  
  /** @brief Function that return beta value */
  Double_t Beta() const { return fP / fE; }

  /** @brief Function that return charge value */
  Int_t Charge() const { return fCharge; }

  /** @brief Function that return curvature */
  Double_t Curvature() const;

  /** @brief Function that return pointer to X direction value */
  const Double_t *Direction() const { return &fXdir; }

  /** @brief Function that return X direction value */
  Double_t DirX() const { return fXdir; }

  /** @brief Function that return Y direction value */
  Double_t DirY() const { return fYdir; }

  /** @brief Function that return Z direction value */
  Double_t DirZ() const { return fZdir; }

  /** @brief Function that return energy value */
  Double_t E() const { return fE; }

  /** @brief Function that return energy deposition value */
  Double_t Edep() const { return fEdep; }

  /** @brief Function that return event number */
  Int_t Event() const { return fEvent; }

  /** @brief Function that return slot number */
  Int_t EventSlot() const { return fEvslot; }

  /** @brief Function that return true if starting from boundary */
  Bool_t FromBoundary() const { return fFrombdr; }

   /** @brief Function that return G5 particle code */
  Int_t G5code() const { return fG5code; }

   /** @brief Function that return element index */
  Int_t EIndex() const { return fEindex; }

  /** @brief Function that return gamma value*/
  Double_t Gamma() const { return fMass ? fE / fMass : TMath::Limits<double>::Max(); }

  /** @brief Function that return selected physical step */
  Double_t GetPstep() const { return fPstep; }

  /** @brief Function that return volume */
  TGeoVolume *GetVolume() const;

  /** @brief Function that return next volume */
  TGeoVolume *GetNextVolume() const;

  /** @brief Function that return material */
  TGeoMaterial *GetMaterial() const;

  /** @brief Function that return path in volume */
  VolumePath_t *GetPath() const { return fPath; }

  /** @brief Function that return next path in volume */
  VolumePath_t *GetNextPath() const { return fNextpath; }
  
  /** @brief Function that return number of physical step made */
  Int_t GetNsteps() const { return fNsteps; }

   /** @brief Function that return physical step */
  Double_t GetStep() const { return fStep; }

  /** @brief Function that return straight distance to next boundary */
  Double_t GetSnext() const { return fSnext; }

  /** @brief Function that return safe distance to any boundary */
  Double_t GetSafety() const { return fSafety; }

  /** @brief Function that check if track is alive */
  Bool_t IsAlive() const { return (fStatus != kKilled); }

   /** @brief Function that check if track is on boundary */
  Bool_t IsOnBoundary() const { return (fStatus == kBoundary); }

  /** @brief Function that return current volume index */
  Int_t Vindex() const { return fVindex; }

  /** @brief Function that set status killed to track */
  void Kill() { fStatus = kKilled; }
  
  /** @brief Function that return mass value */
  Double_t Mass() const { return fMass; }

  /** @brief Function that return momentum value */
  Double_t P() const { return fP; }

  /** @brief Function that return momentum X component */
  Double_t Px() const { return fP * fXdir; }

  /** @brief Function that return momentum Y component */
  Double_t Py() const { return fP * fYdir; }

  /** @brief Function that return momentum Z component */
  Double_t Pz() const { return fP * fZdir; }

  /** @brief Function that return module momentum's value */
  Double_t Pt() const { return fP * Math::Sqrt(fXdir * fXdir + fYdir * fYdir); }
  
  /** @brief Function that return index of corresponding particle */
  Int_t Particle() const { return fParticle; }

  /** @brief Function that set status pending to track */
  Bool_t Pending() const { return fPending; }

  /** @brief Function that return particle pdg code */
  Int_t PDG() const { return fPDG; }

  /** @brief Function that return current process */
  Int_t Process() const { return fProcess; }
  const Double_t *Position() const { return &fXpos; }

  /** @brief Function that return X position */
  Double_t PosX() const { return fXpos; }

  /** @brief Function that return Y position */
  Double_t PosY() const { return fYpos; }

  /** @brief Function that return Z position */
  Double_t PosZ() const { return fZpos; }

  /** @brief Print function */
  void Print(Int_t trackindex = 0) const;

  /** Function that return particle species */
  Species_t Species() const { return fSpecies; }

  /** Function that return track status */
  TrackStatus_t Status() const { return fStatus; }

  /** Function that return time */
  Double_t Time() const { return fTime; }
  
  /** Clear function */
  void Clear(Option_t *option = "");

  /** @brief Function that return X coordinate */
  Double_t X() const { return fXpos; }

  /** @brief Function that return Y coordinate */
  Double_t Y() const { return fYpos; }

  /** @brief Function that return Z coordinate */
  Double_t Z() const { return fZpos; }
  
  /**
   * @brief Function that read track from vector
   * 
   * @param arr Array of tracks
   * @param i Position to read
   */
  void ReadFromVector(const GeantTrack_v &arr, Int_t i);

  /**
   * @brief Function that set event number
   * 
   * @param event Event that should be set as fEvent
   */
  void SetEvent(Int_t event) { fEvent = event; }

  /**
   * @brief Function that set event slot number
   * 
   * @param slot Event slot that should be set as fEvslot
   */
  void SetEvslot(Int_t slot) { fEvslot = slot; }

  /**
   * @brief Function that set particle index
   * 
   * @param particle Particle that should be set as fParticle
   */
  void SetParticle(Int_t particle) { fParticle = particle; }

  /**
   * @brief Function that set particle pdg code
   * 
   * @param pdg Particle pdg code that should be set as fPDG 
   */
  void SetPDG(Int_t pdg) { fPDG = pdg; }

  /**
   * @brief Function that set G5 particle code
   * 
   * @param g5code G5 particle code that should be set as fG5code
   */
  void SetG5code(Int_t g5code) { fG5code = g5code; }

  /**
   * @brief Function that set element index
   * 
   * @param ind Element index that should be set as fEindex
   */
  void SetEindex(Int_t ind) { fEindex = ind; }

  /**
   * @brief Function that set charge
   * 
   * @param charge Charge that should be set as fCharge
   */
  void SetCharge(Int_t charge) { fCharge = charge; }

  /**
   * @brief Function that set process 
   * 
   * @param process Process that should be set as fProcess
   */
  void SetProcess(Int_t process) { fProcess = process; }

  /**
   * @brief Function that set current volume index
   * 
   * @param ind Current volume index that should be set as fVindex
   */
  void SetVindex(Int_t ind) { fVindex = ind; }

  /**
   * @brief Function that set current step
   * 
   * @param nsteps Current step hat should be set as fNsteps
   */
  void SetNsteps(Int_t nsteps) { fNsteps = nsteps; }

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
  void SetMass(Double_t mass) { fMass = mass; }

  /**
   * @brief Function that set X, Y, Z positions
   * 
   * @param x X position
   * @param y Y position
   * @param z Z position
   */
  void SetPosition(Double_t x, Double_t y, Double_t z) {
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
  void SetDirection(Double_t dx, Double_t dy, Double_t dz) {
    fXdir = dx;
    fYdir = dy;
    fZdir = dz;
  }

  /**
   * @brief Function that set momentum
   * 
   * @param p Current momentum should be set as fP
   */
  void SetP(Double_t p) { fP = p; }

  /**
   * @brief Function that set energy
   * 
   * @param e Current E should be set as fE
   */
  void SetE(Double_t e) { fE = e; }

  /**
   * @brief Function that set time
   * 
   * @param time Current time should be set as fTime
   */
  void SetTime(Double_t time) { fTime = time; }

  /**
   * @brief Function that set energy deposition
   * 
   * @param edep Current energy deposition should be set as fEdep
   */
  void SetEdep(Double_t edep) { fEdep = edep; }

  /**
   * @brief Function that set current physical step
   * 
   * @param pstep Current physical step should be set as fPstep
   */
  void SetPstep(Double_t pstep) { fPstep = pstep; }

  /**
   * @brief Function that set straight distance to next boundary
   * 
   * @param snext Straight distance to next boundary should be set as fSnext
   */
  void SetSnext(Double_t snext) { fSnext = snext; }

  /**
   * @brief Function that set safe distance to any boundary
   * 
   * @param safety Safe distance to any boundary hould be set as fSafety
   */
  void SetSafety(Double_t safety) { fSafety = safety; }

  /**
   * @brief Function that set starting from boundary flag
   * 
   * @param flag Flag that is true if starting from boundary 
   */
  void SetFrombdr(Bool_t flag) { fFrombdr = flag; }

  /**
   * @brief Function that set pending status
   * 
   * @param flag Flag that should be set pending
   */
  void SetPending(Bool_t flag) { fPending = flag; }

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

  ClassDefNV(GeantTrack, 1) // The track
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
  cacheline_pad_t pad0_;
  Int_t fMaxtracks;  /** Max size for tracks */
  Int_t fNselected;  /** Number of selected tracks */
  BitSet *fHoles;    /** Bits of holes */
  BitSet *fSelected; /** Mask of selected tracks for the current operation */
  Bool_t fCompact;   /** Flag marking the compactness */
  Bool_t fMixed;     /** Contains tracks in mixed volumes */

#ifdef __STAT_DEBUG_TRK
  GeantTrackStat fStat; /** Statistics for the track container */
#endif
  Int_t fMaxDepth; /** Maximum geometry depth allowed */
  size_t fBufSize; /** Size of the internal buffer */
  char *fVPstart;  /** Address of volume path buffer */
  char *fBuf;      /** Buffer holding tracks data */

  Int_t *fEventV;          /** Event numbers */
  Int_t *fEvslotV;         /** Event slots */
  Int_t *fParticleV;       /** Indices of corresponding particles */
  Int_t *fPDGV;            /** Particle pdg codes */
  Int_t *fG5codeV;         /** G5 internal codes */
  Int_t *fEindexV;         /** Element indices */
  Int_t *fChargeV;         /** Particle charges */
  Int_t *fProcessV;        /** Current process */
  Int_t *fVindexV;         /** Volume index */
  Int_t *fNstepsV;         /** Number of steps made */
  Species_t *fSpeciesV;    /** Particle species */
  TrackStatus_t *fStatusV; /** Track statuses */
  Double_t *fMassV;        /** Particle masses */
  Double_t *fXposV;        /** Array of track X positions */
  Double_t *fYposV;        /** Array of track Y positions */
  Double_t *fZposV;        /** Array of track Z positions */
  Double_t *fXdirV;        /** Array of track directions */
  Double_t *fYdirV;
  Double_t *fZdirV;
  Double_t *fPV;      /** Momentum */
  Double_t *fEV;      /** Energies */
  Double_t *fTimeV;   /** Time */
  Double_t *fEdepV;   /** Energy depositions */
  Double_t *fPstepV;  /** Selected physical steps */
  Double_t *fStepV;   /** Current steps */
  Double_t *fSnextV;  /** Straight distances to next boundary */
  Double_t *fSafetyV; /** Safe distances to any boundary */
  Bool_t *fFrombdrV;  /** True if starting from boundary */
  Bool_t *fPendingV;
  VolumePath_t **fPathV;     /** Paths for the particles in the geometry */
  VolumePath_t **fNextpathV; /** Paths for next volumes */
  
  /**
   * @brief Function that assign in current buffer
   * 
   * @param buff Buffer to be assigned to
   * @param size Size of buffer
   */
  void AssignInBuffer(char *buff, Int_t size);

  /**
   * @brief Function that copy to current buffer
   * 
   * @param buff Buffer to be copied to
   * @param size Size of buffer
   */
  void CopyToBuffer(char *buff, Int_t size);

private:

  /** @brief Copy constructor (not allowed) */
  GeantTrack_v(const GeantTrack_v &track_v);

  /** @brief Operator = */
  GeantTrack_v &operator=(const GeantTrack_v &track_v);

public:

  /** @brief GeantTrack_v constructor */
  GeantTrack_v();

  /** 
   * @brief GeantTrack_v parametrized constructor
   * 
   * @param size Size of track
   * @param maxdepth Maximum allowed geometry depth
   */
  GeantTrack_v(Int_t size, Int_t maxdepth);

  /** @brief GeantTrack_v destructor */
  virtual ~GeantTrack_v();
  
  /** @brief  Function that returned buffer size  */
  size_t BufferSize() const { return fBufSize; }

    /** @brief  Function that returned max size for tracks */
  Int_t Capacity() const { return fMaxtracks; }

  /**
   * @brief Function that compare 2 tracks
   * 
   * @param tr1 First track
   * @param i1 Bit number 'i1'
   * @param tr2 Second track
   * @param i2 Bit number 'i2'
   */
  static Bool_t IsSame(const GeantTrack_v &tr1, Int_t i1, const GeantTrack_v &tr2, Int_t i2);
  
  /**
   * @brief Implementation of memcpy skipping the alignment check.
   */
//   void *memcpy_align(void *dst, const void *src, size_t len) {return memcpy(dst,src,len);}
  static void *memcpy_align(void *dst, const void *src, size_t len) 
  __attribute__((always_inline)) 
  {
    return memcpy(dst,src,len);
    size_t i;
    long *d = (long *)dst;
    const long *s = (const long *)src;
    for (i=0; i<1+len/sizeof(long); ++i)
      d[i] = s[i];
    return dst;
  }    

#ifdef GEANT_CUDA_DEVICE_BUILD
  GEANT_CUDA_BOTH_CODE

  /** @brief  Function that returned number of tracks contained  */
  Int_t GetNtracks() const { return fNtracks; }

  /** @brief  Function that set number of tracks contained  */
  void SetNtracks(Int_t ntracks) { fNtracks = ntracks; }
#else

  /** @brief  Function that returned number of tracks contained  C++11 */
  Int_t GetNtracks() const { return fNtracks.load(); }

   /** @brief  Function that set number of tracks contained  C++11 */
  void SetNtracks(Int_t ntracks) { fNtracks.store(ntracks); }
#endif

   /** @brief  Function that return number of selected tracks  */
  Int_t GetNselected() const { return fNselected; }
#ifdef __STAT_DEBUG_TRK

  /** @brief  Function that return track statistics */
  GeantTrackStat &GetTrackStat() { return fStat; }
#endif
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Add track function
   * 
   * @param track Track that should be added
   * @param import Flag for importing (by default False)
   */
  Int_t AddTrack(GeantTrack &track, Bool_t import = kFALSE);

  /**
   * @brief Add & synctrack function
   * 
   * @param track Track that should be added
   */
  Int_t AddTrackSync(GeantTrack &track);
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Add track function
   * 
   * @param arr Array of tracks 
   * @param i  Bit number 'i'
   * @param import Flag for importing (by default False)
   */
  Int_t AddTrack(GeantTrack_v &arr, Int_t i, Bool_t import = kFALSE);

  /**
   * @brief Add & synctrack function
   * 
   * @param arr Track array
   * @param i Bit number 'i'
   */
  Int_t AddTrackSync(GeantTrack_v &arr, Int_t i);

  /**
   * @brief Add track function
   * 
   * @param arr Track array
   * @param istart Start index of tracks (Start bit number 'i')
   * @param iend End index of tracks (End bit number 'i')
   * @param import Import flag (by default False) 
   */
  void AddTracks(GeantTrack_v &arr, Int_t istart, Int_t iend, Bool_t import = kFALSE);
  
  /** @brief Function that check track */
  void CheckTracks();
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief  Function that mark removed bit number 'i' through value of bits of holes and flag marking the compactness
   * 
   * @param i Bit number 'i'
   */
  void MarkRemoved(Int_t i) {
    fHoles->SetBitNumber(i);
    fCompact = kFALSE;
  }

  /**
   * @brief Function that removes tracks
   * 
   * @param from Start index of tracks to be removed
   * @param to End index of tracks to be removed
   */
  void RemoveTracks(Int_t from, Int_t to);

  /** @brief Function that return size of track */
  size_t Sizeof() const { return sizeof(GeantTrack_v) + fBufSize; }

  /**
   * @brief Function that delete track
   * 
   * @param itr Track ID
   */
  void DeleteTrack(Int_t itr);

  /**
   * @brief Function that deselect bit number 'i' from value
   * 
   * @param i Bit number i
   */
  void Deselect(Int_t i) { fSelected->SetBitNumber(i, kFALSE); }

  /** @brief Function that deselect all bit number from value  */
  void DeselectAll() {
    fSelected->ResetAllBits();
    fNselected = 0;
  }

  /**
   * @brief Function that select bit number 'i' to be value
   * 
   * @param i Bit number i
   */
  void Select(Int_t i) { fSelected->SetBitNumber(i); }

  /**
   * @brief  Function that select tracks
   * @details Function that select tracks
   * 
   * @param n Number of tracks to be selected
   */
  void SelectTracks(Int_t n) { fNselected = n; }

  /**
   * @brief Set if tracks are in mixed volumes through fMixed flag
   * 
   * @param flag Flag that shows if it is contained tracks in mixed volumes
   */
  void SetMixed(Bool_t flag) { fMixed = flag; }

  /**
   * @brief Sorting function for track status
   * 
   * @param status Track status to be checked
   */
  Int_t SortByStatus(TrackStatus_t status);

  /**
   * @brief Sorting function for tracks where the step was limited by discrete processes
   * 
   * @return Number of selected tracks
   */
  Int_t SortByLimitingDiscreteProcess();
  /**
   * @brief Function for removal tracks according status
   * 
   * @param status Track status to be selected for removal
   * @param output Output array of tracks 
   */
  Int_t RemoveByStatus(TrackStatus_t status, GeantTrack_v &output);

  /**
   * @brief Selection function 
   * 
   * @param i Bit number 'i' to be value
   */
  Bool_t IsSelected(Int_t i) { return fSelected->TestBitNumber(i); }

  /** @brief Clear function */
  void Clear(Option_t *option = "");

  /** 
   * @brief Compact function
   * @param  moveto ???????
   */
  Int_t Compact(GeantTrack_v *moveto = 0);

  /**
   * @brief Contain function
   * 
   * @param evstart Event to start from
   * @param nevents Numbe rof evets (by default 1)
   */
  Bool_t Contains(Int_t evstart, Int_t nevents = 1) const;

  /** @brief Clear selection of bit value */
  void ClearSelection() { fSelected->ResetAllBits(); }

  /**
   * @brief Function that return track
   * 
   * @param i Bit number 'i'
   * @param track Track to be returned
   */
  void GetTrack(Int_t i, GeantTrack &track) const;

  /** @brief Function that check flag marking the compactness */
  Bool_t IsCompact() const { return fCompact; }

  /** @brief Function that return how many tracks in mixed volumes */
  Bool_t IsMixed() const { return fMixed; }

  /** @brief Function that print pointers */
  void PrintPointers() { printf("fEventV=%p fFrombdrV=%p\n", (void *)fEventV, (void *)fFrombdrV); }

  /**
   * @brief Function that print track
   * 
   * @param itr Track ID
   */
  void PrintTrack(Int_t itr) const;

  /** @brief Function that print all tracks */
  void PrintTracks() const;

  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function for navigation that find next boundary and step
   * 
   * @param ntracks Number of tracks
   * @param pstep Previos step
   * @param x X position
   * @param y Y position
   * @param z Z position
   * @param dirx X direction
   * @param diry Y direction
   * @param dirz Z direction
   * @param pathin Path inside in the volume
   * @param pathout Path outside in the volume
   * @param step Step to be proccessed
   * @param safe Safety distance
   * @param isonbdr 
   * @param trk Track
   */
  void NavFindNextBoundaryAndStep(Int_t ntracks, const Double_t *pstep, const Double_t *x,
                                  const Double_t *y, const Double_t *z, const Double_t *dirx,
                                  const Double_t *diry, const Double_t *dirz, VolumePath_t **pathin,
                                  VolumePath_t **pathout, Double_t *step, Double_t *safe,
                                  Bool_t *isonbdr, const GeantTrack_v *trk);
  
  /**
   * @brief Function for navigation that check if location is the same or not
   * 
   * @param ntracks Number of tracks
   * @param start Start volume path
   * @param end End volume path
   * @param same Boolean flag that check same location
   */
  void NavIsSameLocation(Int_t ntracks, VolumePath_t **start, VolumePath_t **end, Bool_t *same);
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function for navigation that check if location is the same or not for single track
   * 
   * @param itr Track ID
   * @param start Start volume path
   * @param end End volume path
   */
  Bool_t NavIsSameLocationSingle(Int_t itr, VolumePath_t **start, VolumePath_t **end);

// void InspectGeometryState(Int_t itr) const;
// void InspectIsSameLocation(Int_t itr) const;

#ifdef USE_VECGEOM_NAVIGATOR

  /**
   * @brief Function that check if location path's consistency
   * 
   * @param itr Track ID
   */
  void CheckLocationPathConsistency(Int_t itr) const;
#endif

  /**
   * @brief Function that provides postponed action for tracks
   * 
   * @param ntracks Number of tracks
   */
  TransportAction_t PostponedAction(Int_t ntracks) const;
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function that provides postponed action for track
   * 
   * @param itr Track ID
   * @param output Output of tracks 
   */
  Int_t PostponeTrack(Int_t itr, GeantTrack_v &output);

  /**
   * @brief Function that provides postponed action for all track
   * 
   * @param output Output of tracks
   */
  Int_t PostponeTracks(GeantTrack_v &output);
  // void      PropagateBack(Int_t itr, Double_t crtstep);
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function that compute transport length
   * 
   * @param ntracks Number of tracks
   */
  void ComputeTransportLength(Int_t ntracks);
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function that compute single transport length ?????
   * 
   * @param itr Track ID
   */
  void ComputeTransportLengthSingle(Int_t itr);

  /**
   * @brief Function of propagation in volume
   * 
   * @param ntracks Number of tracks
   * @param crtstep ??????
   * @param tid Track ID 
   */
  void PropagateInVolume(Int_t ntracks, const Double_t *crtstep, Int_t tid);
  GEANT_CUDA_BOTH_CODE

  /**
   * @brief Function of propagation of track in volume
   * 
   * @param i Bit number 'i'
   * @param crtstep ???????
   * @param tid Track ID
   */
  void PropagateInVolumeSingle(Int_t i, Double_t crtstep, Int_t tid);

  /**
   * @brief Popagation function in straight trajectories 
   * 
   * @param ntracks Number of tracks
   * @param crtstep  ?????
   */
  Int_t PropagateStraight(Int_t ntracks, Double_t *crtstep);

  /**
   * @brief Function of propagation of tracks
   * 
   * @param output Output array of tracks
   * @param tid Track ID
   */
  Int_t PropagateTracks(GeantTrack_v &output, Int_t tid);
  GEANT_CUDA_BOTH_CODE


  Int_t PropagateTracksSingle(GeantTrack_v &output, Int_t tid, Int_t stage = 0);

  /** 
   * @brief Resize function
   * @param newsize New size to be resized
   */
  void Resize(Int_t newsize);

  /** @brief Function that replace track at positions i with track on position j */
  void ReplaceTrack(Int_t i, Int_t withj);

  /** @brief Reshuffle function */
  Int_t Reshuffle();

  /** 
   * @brief Function that swap tracks at positions i and j
   * @param i Input bit number 'i'
   * @param j Input bit number 'j'
   */
  void SwapTracks(Int_t i, Int_t j);
  
  /** 
   * @brief Function that return beta value
   * @param  i Input bit number 'i'
   */
  Double_t Beta(Int_t i) const { return fPV[i] / fEV[i]; }
  GEANT_CUDA_BOTH_CODE

  /** 
   * @brief Function that return curvature in different areas of geometry
   * @param  i Input bit number 'i'
   */
  Double_t Curvature(Int_t i) const;
  GEANT_CUDA_BOTH_CODE

  /** @brief Function that return safe length */
  Double_t SafeLength(Int_t i, Double_t eps = 1.E-4);

  /** 
   * @brief Function that return gamma value
   * @param  i Input bit number 'i'
   */
  Double_t Gamma(Int_t i) const {
    return fMassV[i] ? fEV[i] / fMassV[i] : TMath::Limits<double>::Max();
  }

  /** 
   * @brief Function that return X projection of momentum value
   * @param  i Input bit number 'i'
   */
  Double_t Px(Int_t i) const { return fPV[i] * fXdirV[i]; }

  /** 
   * @brief Function that return Y projection of momentum value
   * @param  i Input bit number 'i'
   */
  Double_t Py(Int_t i) const { return fPV[i] * fYdirV[i]; }
  GEANT_CUDA_BOTH_CODE

  /** 
   * @brief Function that return Z projection of momentum value
   * @param  i Input bit number 'i'
   */
  Double_t Pz(Int_t i) const { return fPV[i] * fZdirV[i]; }
  GEANT_CUDA_BOTH_CODE

  /** 
   * @brief Function that return module of momentum value 
   * @param  i Input bit number 'i'
   */
  Double_t Pt(Int_t i) const {
    return fPV[i] * Math::Sqrt(fXdirV[i] * fXdirV[i] + fYdirV[i] * fYdirV[i]);
  }

 /** 
  * @brief Function that returnes TGeoVolume
  * @param  i Input bit number 'i'
  */
  TGeoVolume *GetVolume(Int_t i) const;

 /** 
  * @brief Function that returnes next TGeoVolume
  * @param  i Input bit number 'i'
  */
  TGeoVolume *GetNextVolume(Int_t i) const;

  /** 
   * @brief Function that returnes TGeoMaterial 
   * @param  i Input bit number 'i'
   */
  TGeoMaterial *GetMaterial(Int_t i) const;

  /**
   * @brief Function round up align ?????
   * @param num Number ?????
   */
  static Int_t round_up_align(Int_t num) {
    int remainder = num % ALIGN_PADDING;
    if (remainder == 0)
      return num;
    return (num + ALIGN_PADDING - remainder);
  }

  ClassDefNV(GeantTrack_v, 1) // SOA for GeantTrack class
};

#endif
