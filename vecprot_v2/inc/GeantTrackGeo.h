//===--- GeantTrackGeo.h - GeantV ---------------------------------*- C++ -*-===//
//
//                     GeantV Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTrackGeo.h
 * @brief Non-concurrent SOA for geometry tracks.
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TRACK_GEO
#define GEANT_TRACK_GEO

#include <vector>
#include "GeantTrackVec.h"

namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {


/**
 * @brief AOS for tracs used at processing time by geometry
 *
 */
struct GeantTrackGeo {
  GeantTrack *fOriginal; /** Original track from which this was extracted */
  double fXpos;          /** X position */
  double fYpos;          /** Y position */
  double fZpos;          /** Z position */
  double fXdir;          /** X direction */
  double fYdir;          /** Y direction */
  double fZdir;          /** Z direction */
  double fPstep;         /** Selected physical step */
  double fStep;          /** Current step */
  double fSnext;         /** Straight distance to next boundary */
  double fSafety;        /** Safe distance to any boundary */
  bool   fBoundary;      /** True if starting from boundary */
  VolumePath_t *fPath;   /** Path to current volume containing the particle */
  VolumePath_t *fNextpath; /** Path to volume particle is entering into */
};

/**
 * @brief SOA for tracs used at processing time by geometry
 *
 */
class GeantTrackGeo_v {
  using TrackArray_t = std::vector<GeantTrack*>;

public:
  int fNtracks;      /** Number of tracks contained */
  int fMaxtracks;    /** Max size for tracks */

  size_t fBufSize; /** Size of the internal buffer */
  char *fBuf;      /** Buffer holding tracks data */

  GeantTrack **fOriginalV; /** Track originals */
  double *fXposV;          /** Array of track X positions */
  double *fYposV;          /** Array of track Y positions */
  double *fZposV;          /** Array of track Z positions */
  double *fXdirV;          /** Array of track directions */
  double *fYdirV;
  double *fZdirV;
  double *fPstepV;  /** Selected physical steps */
  double *fStepV;   /** Current steps */
  double *fSnextV;  /** Straight distances to next boundary */
  double *fSafetyV; /** Safe distances to any boundary */
  bool *fBoundaryV; /** True if starting from boundary */
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

  GeantTrackGeo_v(const GeantTrackGeo_v &track_v) = delete;
  GeantTrackGeo_v &operator=(const GeantTrackGeo_v &track_v) = delete;

  /**
   * @brief GeantTrackGeo_v constructor based on a provided single buffer.
   */
  GEANT_CUDA_BOTH_CODE
  GeantTrackGeo_v(void *addr, unsigned int nTracks);

public:
  /** @brief GeantTrackGeo_v default constructor */
  GeantTrackGeo_v();

  /**
   * @brief GeantTrackGeo_v parametrized constructor
   *
   * @param size Initial capacity
   */
  GeantTrackGeo_v(int size);

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  GEANT_CUDA_BOTH_CODE
  static GeantTrackGeo_v *MakeInstanceAt(void *addr, unsigned int nTracks);

  /** @brief GeantTrackGeo_v destructor */
  ~GeantTrackGeo_v();

  /** @brief return the contiguous memory size needed to hold a GeantTrackGeo_v */
  GEANT_CUDA_BOTH_CODE
  static size_t SizeOfInstance(size_t nTracks);

  /** @brief  Function that returns the buffer size  */
  size_t BufferSize() const { return fBufSize; }

  /** @brief  Return the address of the data memory buffer  */
  void *Buffer() const { return fBuf; }

  /** @brief  Function that returns buffer size needed to hold the data for nTracks*/
  GEANT_CUDA_BOTH_CODE
  static size_t BufferSize(size_t nTracks);

  /** @brief  Function that returned max size for tracks */
  int Capacity() const { return fMaxtracks; }

  /** @brief  Check direction normalization within tolerance */
  bool IsNormalized(int itr, double tolerance = 1.E-8) const;

  /** @brief  Function that returned number of tracks contained  */
  GEANT_CUDA_BOTH_CODE
  int GetNtracks() const { return fNtracks; }

  /**
   * @brief Add track function
   *
   * @param track Track that should be added
   * @param import Flag for importing (by default False)
   */
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  int AddTrack(GeantTrack &track) {
    int itrack = fNtracks;
    if (itrack == fMaxtracks) {
#ifndef GEANT_CUDA_DEVICE_BUILD
      Resize(2 * fMaxtracks);
#endif
      printf("Error in GeantTrackGeo::AddTrack, resizing is not supported in device code\n");
#endif
    }
    fOriginalV[itrack] = &track;
    fXposV[itrack] = track.fXpos;
    fYposV[itrack] = track.fYpos;
    fZposV[itrack] = track.fZpos;
    fXdirV[itrack] = track.fXdir;
    fYdirV[itrack] = track.fYdir;
    fZdirV[itrack] = track.fZdir;
    fPstepV[itrack] = track.fPstep;
    fStepV[itrack] = track.fStep;
    fSnextV[itrack] = track.fSnext;
    fSafetyV[itrack] = track.fSafety;
    fBoundaryV[itrack] = track.fBoundary;
    fPathV[itrack] = track.fPath;
    fNextpathV[itrack] = track.fNextpath;
    fNtracks++;
    return itrack;
  }
  
  /**
   * @brief Add tracks from a AOS vector into the SOA
   *
   * @param array Array of tracks that should be added
   */
  GEANT_CUDA_BOTH_CODE
  int AddTracks(TrackArray_t const &array);

  /**
   * @brief Update a single original track from the container
   *
   * @param itr Track to update
  */
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  void UpdateOriginalTrack(int itr) const {
    // Update all the original tracks. This should ideally vectorize.
    for (int itr=0; itr<fNtracks; ++itr) {
      GeantTrackGeo &track = *fOriginals[itr];
      track->fXpos = fXposV[itr];
      track->fYpos = fYposV[itr];
      track->fZpos = fZposV[itr];
      track->fXdir = fXdirV[itr];
      track->fYdir = fYdirV[itr];
      track->fZdir = fZdirV[itr];
      track->fPstep = fPstepV[itr];
      track->fStep = fStepV[itr];
      track->fSnext = fSnextV[itr];
      track->fSafety = fSafetyV[itr];
      track->fBoundary = fBoundaryV[itr];
      // The path and nextpath members are already pointing to the storage used by the original track
    }
  }

  /**
   * @brief Update all original tracks from the container
  */
  GEANT_CUDA_BOTH_CODE
  void UpdateOriginalTracks() const;

  /** @brief Function that return size of track */
  size_t Sizeof() const { return sizeof(GeantTrackGeo_v) + fBufSize; }

  /** @brief Function to normalize direction */
  void Normalize(int itr) __attribute__((always_inline)) {
    double norm = 1. / Math::Sqrt(fXdirV[itr] * fXdirV[itr] + fYdirV[itr] * fYdirV[itr] + fZdirV[itr] * fZdirV[itr]);
    fXdirV[itr] *= norm;
    fYdirV[itr] *= norm;
    fZdirV[itr] *= norm;
  }

  /** @brief Clear function */
  GEANT_CUDA_BOTH_CODE
  GEANT_INLINE
  void Clear() { fNtracks = 0; }

  /**
   * @brief Function that print track
   *
   * @param itr Track ID
   */
  GEANT_CUDA_BOTH_CODE
  void PrintTrack(int itr, const char *msg = "") const;

  /** @brief Function that print all tracks */
  void PrintTracks(const char *msg = "") const;

#ifdef USE_VECGEOM_NAVIGATOR

  /**
   * @brief Function that check if location path's consistency
   *
   * @param itr Track ID
   */
  void CheckLocationPathConsistency(int itr) const;
#endif

  /**
   * @brief Check consistency of track navigation
   * @param  itr Track number to be checked
   */
  bool CheckNavConsistency(int itr);

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

  //  ClassDefNV(GeantTrackGeo_v, 1) // SOA for GeantTrack class
};
} // GEANT_IMPL_NAMESPACE

#ifdef GEANT_CUDA
#ifdef GEANT_NVCC
namespace cxx {
class GeantTrackGeo_v;
}
#else
namespace cuda {
class GeantTrackGeo_v;
}
#endif

bool ToDevice(vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo_v> dest, cxx::GeantTrackGeo_v *source, cudaStream_t stream);
bool FromDevice(cxx::GeantTrackGeo_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo_v> source, cudaStream_t stream);
void FromDeviceConversion(cxx::GeantTrackGeo_v *dest, vecgeom::cxx::DevicePtr<cuda::GeantTrackGeo_v> source);
#endif

} // Geant

#endif
