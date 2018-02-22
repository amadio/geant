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
#include <Track.h>

namespace geant {
inline namespace GEANT_IMPL_NAMESPACE {

/**
 * @brief AOS for tracs used at processing time by geometry
 *
 */
struct GeantTrackGeo {
  Track *fOriginal; /** Original track from which this was extracted */
  size_t fId;            /** Id of original track in its container */
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
  bool   fCompSafety;    /** Safety needs to be computed (not on boundary) */
};

/**
 * @brief SOA for tracs used at processing time by geometry
 *
 */
class GeantTrackGeo_v {
#ifndef VECCORE_CUDA
  typedef std::vector<Track *> TrackVec_t;
#else
  typedef vecgeom::Vector<Track *> TrackVec_t;
#endif

public:
  int fNtracks;      /** Number of tracks contained */
  int fMaxtracks;    /** Max size for tracks */

  size_t fBufSize; /** Size of the internal buffer */
  char *fBuf;      /** Buffer holding tracks data */

  Track **fOriginalV; /** Track originals */
  size_t *fIdV;            /** Id's in the original container */
  double *fXposV;          /** Arrays of track positions */
  double *fYposV;
  double *fZposV;
  double *fXdirV;          /** Arrays of track directions */
  double *fYdirV;
  double *fZdirV;
  double *fPstepV;         /** Selected physical steps */
  double *fStepV;          /** Current steps */
  double *fSnextV;         /** Straight distances to next boundary */
  double *fSafetyV;        /** Safe distances to any boundary */
  bool   *fCompSafetyV;    /** Safety needs to be computed (not on boundary) */

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

  GeantTrackGeo_v(const GeantTrackGeo_v &track_v) = delete;
  GeantTrackGeo_v &operator=(const GeantTrackGeo_v &track_v) = delete;

  /**
   * @brief GeantTrackGeo_v constructor based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
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
   * @brief Track MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_DEVICE
  static GeantTrackGeo_v *MakeInstanceAt(void *addr, unsigned int nTracks);

  /** @brief GeantTrackGeo_v destructor */
  ~GeantTrackGeo_v();

  /** @brief return the contiguous memory size needed to hold a GeantTrackGeo_v */
  VECCORE_ATT_DEVICE
  static size_t SizeOfInstance(size_t nTracks);

  /** @brief  Function that returns the buffer size  */
  size_t BufferSize() const { return fBufSize; }

  /** @brief  Return the address of the data memory buffer  */
  void *Buffer() const { return fBuf; }

  /** @brief  Function that returns buffer size needed to hold the data for nTracks*/
  VECCORE_ATT_HOST_DEVICE
  static size_t BufferSize(size_t nTracks);

  /** @brief  Function that returned max size for tracks */
  int Capacity() const { return fMaxtracks; }

  /** @brief  Check direction normalization within tolerance */
  bool IsNormalized(int itr, double tolerance = 1.E-8) const;

  /** @brief  Function that returned number of tracks contained  */
  VECCORE_ATT_HOST_DEVICE
  int GetNtracks() const { return fNtracks; }

  /**
   * @brief Add track function
   *
   * @param track Track that should be added
   * @param import Flag for importing (by default False)
   */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  int AddTrack(Track &track, size_t id) {
    int itrack = fNtracks;
    if (itrack == fMaxtracks) {
#ifndef VECCORE_CUDA_DEVICE_COMPILATION
      Resize(2 * fMaxtracks);
#else
      printf("Error in GeantTrackGeo::AddTrack, resizing is not supported in device code\n");
#endif
    }
    fOriginalV[itrack] = &track;
    fIdV[itrack]   = id;
    fXposV[itrack] = track.X();
    fYposV[itrack] = track.Y();
    fZposV[itrack] = track.Z();
    fXdirV[itrack] = track.Dx();
    fYdirV[itrack] = track.Dy();
    fZdirV[itrack] = track.Dz();
    fPstepV[itrack] = track.GetPstep();
    fStepV[itrack] = track.GetStep();
    fSnextV[itrack] = track.GetSnext();
    fSafetyV[itrack] = track.GetSafety();
    fCompSafetyV[itrack] = !track.Boundary();
    fNtracks++;
    return itrack;
  }
  
  /**
   * @brief Add tracks from a AOS vector into the SOA
   *
   * @param array Array of tracks that should be added
   */
  VECCORE_ATT_DEVICE
  int AddTracks(TrackVec_t const &array);

  /**
   * @brief Update a single original track from the container
   *
   * @param itr Track to update
  */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdateOriginalTrack(int itr) const {
    // Update the original track itr.
    Track &track = *fOriginalV[itr];
    track.SetPosition(fXposV[itr], fYposV[itr], fZposV[itr]);
    track.SetDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
    track.SetPstep(fPstepV[itr]);
    track.SetStep(fStepV[itr]);
    track.SetSnext(fSnextV[itr]);
    track.SetSafety(fSafetyV[itr]);
    track.SetBoundary(!fCompSafetyV[itr]);
  }

  /**
   * @brief Update all original tracks from the container
  */
  VECCORE_ATT_HOST_DEVICE
  void UpdateOriginalTracks() const {
    // Update all the original tracks. This should ideally vectorize.
    for (int itr=0; itr<fNtracks; ++itr) {
      UpdateOriginalTrack(itr);
    }
  }

  /** @brief Update tracks positions */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdatePositions() const {
    for (int itr=0; itr<fNtracks; ++itr)
      fOriginalV[itr]->SetPosition(fXposV[itr], fYposV[itr], fZposV[itr]);
  }

  /** @brief Update tracks directions */
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void UpdateDirections() const {
    for (int itr=0; itr<fNtracks; ++itr)
      fOriginalV[itr]->SetDirection(fXdirV[itr], fYdirV[itr], fZdirV[itr]);
  }
      

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
  VECCORE_ATT_HOST_DEVICE
  GEANT_FORCE_INLINE
  void Clear() { fNtracks = 0; }

  /**
   * @brief Resize function
   * @param newsize New size to be resized
   */
  void Resize(int newsize);

  /**
   * @brief Function rounding up an integer to the aligned value
   * @param num Value be aligned
   */
  VECCORE_ATT_HOST_DEVICE
  static int RoundUpAlign(int num) {
    int remainder = num % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return num;
    return (num + GEANT_ALIGN_PADDING - remainder);
  }

  /**
   * @brief Function rounding up an address to the aligned value
   * @param buf Address to be aligned
   */
  VECCORE_ATT_HOST_DEVICE
  static char *RoundUpAlign(char *buf) {
    long remainder = ((long)buf) % GEANT_ALIGN_PADDING;
    if (remainder == 0)
      return buf;
    return (buf + GEANT_ALIGN_PADDING - remainder);
  }

  //  ClassDefNV(GeantTrackGeo_v, 1) // SOA for Track class
};
} // GEANT_IMPL_NAMESPACE

#ifdef VECCORE_CUDA
namespace cxx {
class GeantTrackGeo_v;
}
#else
namespace cuda {
class GeantTrackGeo_v;
}
#endif

} // Geant

#endif
