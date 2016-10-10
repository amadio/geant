//===--- GeantTaskData.h - Geant-V ----------------------------*- C++ -*-===//
//
//                     Geant-V Prototype
//
//===----------------------------------------------------------------------===//
/**
 * @file GeantTaskData.h
 * @brief Implementation of data organized per thread Geant-V prototype
 * @author Andrei Gheata
 */
//===----------------------------------------------------------------------===//

#ifndef GEANT_TASKDATA
#define GEANT_TASKDATA

#ifndef GEANT_TRACK
#include "GeantTrackVec.h"
#endif
 #include "GeantPropagator.h"

#include <deque>
#include <vector>

#include "Geant/Typedefs.h"

namespace geantphysics {
  class PhysicsData;
}

#ifdef USE_ROOT
class TRandom;
#endif
#ifdef USE_VECGEOM_NAVIGATOR
#include "base/RNG.h"
using VECGEOM_NAMESPACE::RNG;
#endif

#ifdef VECCORE_CUDA
#include "base/Vector.h"
#else
#include <vector>
#endif

/**
 * @brief Class GeantTaskData
 * @details Class descripting data organized per thread
 *
 */
namespace Geant {
inline namespace GEANT_IMPL_NAMESPACE {

class GeantBasketMgr;
class GeantBasket;
class GeantTrackGeo_v;

class GeantTaskData {
public:
#ifdef VECCORE_CUDA
  template <class T>
  using vector_t = vecgeom::Vector<T>;
#else
  template <class T>
  using vector_t = std::vector<T>;
#endif

  GeantPropagator *fPropagator; /** GeantPropagator */
  int fTid;              /** Thread unique id */
  int fNode;             /** Locality node */
  size_t fNthreads;      /** Number of transport threads */
  int fMaxDepth;         /** Maximum geometry depth */
  int fSizeBool;         /** Size of bool array */
  int fSizeDbl;          /** Size of dbl array */
  bool fToClean;         /** Flag set when the basket queue is to be cleaned */
  Volume_t *fVolume;     /** Current volume per thread */
#ifdef USE_VECGEOM_NAVIGATOR
  RNG *fRndm;            /** Random generator for thread */
#elif USE_ROOT
  TRandom *fRndm;        /** Random generator for thread */
#endif
  bool *fBoolArray;      /** [fSizeBool] Thread array of bools */
  double *fDblArray;     /** [fSizeDbl] Thread array of doubles */
  GeantTrack fTrack;     /** Track support for this thread */
  VolumePath_t *fPath;   /** Volume path for the thread */
  VolumePath_t **fPathV;    /** Volume path for the thread */
  VolumePath_t **fNextpathV; /** Volume path for the thread */
  GeantTrackGeo_v *fGeoTrack; /** Geometry track SOA */
  GeantBasketMgr *fBmgr; /** Basket manager collecting mixed tracks */
  GeantBasket *fReused;  /** Basket having tracks to be reused in the same volume */
  GeantBasket *fImported;/** Basket used to import tracks from the event server */
#ifdef VECCORE_CUDA
  char fPool[sizeof(std::deque<GeantBasket *>)]; // Use the same space ...
#else
  std::deque<GeantBasket *> fPool; /** Pool of empty baskets */
#endif
  int fSizeInt;                             // current size of IntArray
  int *fIntArray;                           // Thread array of ints (used in vector navigation)
  GeantTrack_v  *fTransported;              // Transported tracks in current step
  vector_t<GeantTrack *> fTransported1;     // Transported tracks in current step
  int            fNkeepvol;                 // Number of tracks keeping the same volume
  int fNsteps;           /** Total number of steps per thread */
  int fNsnext;           /** Total number of calls to getting distance to next boundary */
  int fNphys;            /** Total number of steps to physics processes */
  int fNmag;             /** Total number of partial steps in magnetic field */
  int fNpart;            /** Total number of particles transported by the thread */
  int fNsmall;           /** Total number of small steps taken */
  int fNcross;           /** Total number of boundary crossings */

  geantphysics::PhysicsData  *fPhysicsData; /** Physics data per thread */

private:
   // a helper function checking internal arrays and allocating more space if necessary
  template <typename T> static void CheckSizeAndAlloc(T *&array, int &currentsize, size_t wantedsize) {
     if (wantedsize < (size_t) currentsize)
      return;
    T *newarray = new T[wantedsize];
    memcpy(newarray,array,currentsize*sizeof(T));
    delete[] array;
    array = newarray;
    currentsize = wantedsize;
  }

  /**
   * @brief GeantTaskData constructor based on a provided single buffer.
   */
  VECCORE_ATT_HOST_DEVICE
  GeantTaskData(void *addr, size_t nTracks, int maxdepth, int maxPerBasket);

public:
  /** @brief GeantTaskData constructor */
  GeantTaskData(size_t nthreads, int maxDepth, int maxPerBasket);

  /** @brief GeantTaskData destructor */
  ~GeantTaskData();

  /**
   * @brief GeantTrack MakeInstance based on a provided single buffer.
   */
  VECCORE_ATT_HOST_DEVICE
  static GeantTaskData *MakeInstanceAt(void *addr, size_t nTracks, int maxdepth, int maxPerBasket);

  /** @brief return the contiguous memory size needed to hold a GeantTrack_v size_t nTracks, size_t maxdepth */
  VECCORE_ATT_HOST_DEVICE
  static size_t SizeOfInstance(size_t nthreads, int maxDepth, int maxPerBasket);

  /**
   * @brief Function that return double array
   *
   * @param size Size of double array
   */
  double *GetDblArray(int size) {
    CheckSizeAndAlloc<double>(fDblArray, fSizeDbl, size);
    return fDblArray;
  }

  /**
   * @brief Function that return boolean array
   *
   * @param size Size of boolean array
   */
  bool *GetBoolArray(int size) {
    CheckSizeAndAlloc<bool>(fBoolArray, fSizeBool, size);
    return fBoolArray;
  }

  /**
   * @brief Function that returns int array
   *
   * @param size Size of int array
   */
  int *GetIntArray(int size) {
    CheckSizeAndAlloc<int>(fIntArray, fSizeInt, size);
    return fIntArray;
  }

  /**
   * @brief Function that returns a (per thread/task) preallocated NavigationState object
   *
   */
  VECCORE_ATT_HOST_DEVICE
  VolumePath_t *GetPath() {
    return fPath;
  }

  /**
   * @brief Get the cleared storedtrack
   */
  GeantTrack &GetTrack() {
    fTrack.Clear();
    return fTrack;
  }

  /**
   * @brief Get next free basket or null if not available
   * @details Get pointer to next free basket
   */
  GeantBasket *GetNextBasket();

  /*
   * @brief Return the size of the basket pool
   *
   */
#ifndef VECCORE_CUDA
  size_t GetBasketPoolSize() const { return fPool.size(); }
#endif

  /** @brief Setter for the toclean flag */
  void SetToClean(bool flag) { fToClean = flag; }

  /** @brief Getter for the toclean flag */
  bool NeedsToClean() const { return fToClean; }

  /**
   * @brief Recycles a given basket
   *
   * @param b Pointer to current GeantBasket for recycling
   */
  void RecycleBasket(GeantBasket *b);

  /**
   * @brief Function cleaning a number of free baskets
   *
   * @param ntoclean Number of baskets to be cleaned
   * @return Number of baskets actually cleaned
   */
  int CleanBaskets(size_t ntoclean);

private:
  /**
   * @brief Constructor GeantTaskData
   * @todo Still not implemented
   */
  GeantTaskData(const GeantTaskData &);

  /**
   * @brief Operator &operator=
   * @todo Still not implemented
   */
  GeantTaskData &operator=(const GeantTaskData &);

  // ClassDef(GeantTaskData, 1) // Stateful data organized per thread
};
} // GEANT_IMPL_NAMESPACE
} // Geant
#endif
