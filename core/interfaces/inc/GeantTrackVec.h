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

#ifndef GEANT_TRACK_VEC
#define GEANT_TRACK_VEC

#include "base/BitSet.h"
#include "GeantTrack.h"
#include "FieldLookup.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantTrack_v;
class GeantTaskData;
class GeantPropagator;


/**
 * @brief SOA for GeantTrack used at processing time
 *
 */
class GeantTrack_v {
public:
  using BitSet = veccore::BitSet;
  static size_t const cacheline_size = 64;
  typedef char cacheline_pad_t[cacheline_size];
#ifdef VECCORE_CUDA_DEVICE_COMPILATION
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
  int *fMotherV;           /** Indices of mother particles */
  int *fPDGV;              /** Particle pdg codes */
  int *fGVcodeV;           /** GV internal codes */
  int *fEindexV;           /** Element indices */
  int *fBindexV;           /** Indices for the track blocks */
  int *fChargeV;           /** Particle charges */
  int *fProcessV;          /** Current process */
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
  static GeantTrack_v *MakeInstanceAt(void *addr, unsigned int nTracks, int maxdepth);

  /** @brief GeantTrack_v destructor */
  ~GeantTrack_v();

  /** @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth */
  VECCORE_ATT_HOST_DEVICE
  static size_t SizeOfInstance(size_t nTracks, size_t maxdepth);

  /** @brief  Function that returns the buffer size  */
  size_t BufferSize() const { return fBufSize; }

  /** @brief  Return the address of the data memory buffer  */
  void *Buffer() const { return fBuf; }

  /** @brief  Function that returns buffer size needed to hold the data for nTracks and maxdepth */
  VECCORE_ATT_HOST_DEVICE
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

#ifdef VECCORE_CUDA_DEVICE_COMPILATION

  /** @brief  Function that returned number of tracks contained  */
  VECCORE_ATT_HOST_DEVICE
  int GetNtracks() const { return fNtracks; }

  /** @brief  Function that set number of tracks contained  */
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
  void ComputeTransportLength(int ntracks, GeantTaskData *);

  /**
   * @brief Function that compute single transport length ?????
   *
   * @param itr Track ID
   */
  VECCORE_ATT_HOST_DEVICE
  void ComputeTransportLengthSingle(int itr, GeantTaskData *);

  /**
   * @brief Function of propagation in volume
   *
   * @param ntracks Number of tracks
   * @param crtstep ??????
   * @param tid Track ID
   */
  VECCORE_ATT_HOST_DEVICE
  void PropagateInVolume(int ntracks, const double *crtstep, GeantTaskData *td);

  /**
   * @brief Function of propagation of track in volume
   *
   * @param i Bit number 'i'
   * @param crtstep ???????
   * @param tid Track ID
   */
  VECCORE_ATT_HOST_DEVICE
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

  VECCORE_ATT_HOST_DEVICE
  int PropagateTracksScalar(GeantTaskData *td, int stage = 0);

  VECCORE_ATT_HOST_DEVICE
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
   * @brief Function that return magnetic field for particle 'i'
   * @param  i             Input track number 'i'
   * @param  BfieldOut[3]  Output magnetic field vector value (global coordinates)
   * @param  bmagOut       Output (optional) field magnitude
   */
  VECCORE_ATT_HOST_DEVICE
  void GetFieldValue( /*GeantTaskData *td,*/ int i, double BfieldOut[3], double *bmagOut) const;
  //  {
  //     vecgeom::Vector3D<double> Position (fXposV[i], fYposV[i], fZposV[i]);
  //     FieldLookup::GetFieldValue(td, Position, BfieldOut, bmagOut);
  //  }

  /** @brief Simple check of tracks -- if a problem is found, then the track is printed 
   * @param  itr           Input track number 'i'
   * @param  msg           Message to inform about location from call is made
   * @param  epsilon       tolerance for problems
   */ 
  VECCORE_ATT_HOST_DEVICE
  void CheckTrack(int itr, const char *msg, double epsilon= 0.0) const;

  /** @brief Check the direction -- report if its length is within epsilon of 1.0
   * @param  itr           Input track number 'i'
   * @param  epsilon       tolerance for problems
   */    
  VECCORE_ATT_HOST_DEVICE
  bool CheckDirection(int itr, double epsilon = 1.0e-6 ) const;
   
  /**
   * @brief Function that return curvature in different areas of geometry
   * @param  i Input bit number 'i'
   */
  VECCORE_ATT_HOST_DEVICE
  double Curvature(int i) const;
   
  /** @brief Function that return safe length */
  VECCORE_ATT_HOST_DEVICE
  double SafeLength(int i, double eps = 1.E-4);

  /**
   * @brief Function that return gamma value
   * @param  i Input bit number 'i'
   */
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
  VECCORE_ATT_HOST_DEVICE
  double Pz(int i) const { return fPV[i] * fZdirV[i]; }

  /**
   * @brief Function that return module of momentum value
   * @param  i Input bit number 'i'
   */
  VECCORE_ATT_HOST_DEVICE
  double Pt(int i) const { return fPV[i] * Math::Sqrt(fXdirV[i] * fXdirV[i] + fYdirV[i] * fYdirV[i]); }

  /** @brief Function that return time traveled in the step */
  VECCORE_ATT_HOST_DEVICE
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
  VECCORE_ATT_HOST_DEVICE
  bool BreakOnStep(int evt, int trk, int stp, int nsteps = 1, const char *msg = "", int itr = -1);

  /**
   * @brief Check consistency of track navigation
   * @param  itr Track number to be checked
   */
  bool CheckNavConsistency(int itr);

  //  ClassDefNV(GeantTrack_v, 1) // SOA for GeantTrack class
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA_ENABLED
#ifdef VECCORE_CUDA
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
